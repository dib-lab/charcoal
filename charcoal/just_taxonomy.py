#! /usr/bin/env python
"""
Remove bad contigs based solely on taxonomy.
"""
import sys
import argparse
import gzip
from collections import Counter, defaultdict
import csv
import os.path

import screed

import sourmash
from sourmash.lca.command_index import load_taxonomy_assignments
from sourmash.lca import LCA_Database, LineagePair

from . import utils
from . import lineage_db
from .lineage_db import LineageDB
from .version import version
from .utils import (get_idents_for_hashval, gather_lca_assignments,
    count_lca_for_assignments, pretty_print_lineage, pretty_print_lineage2,
    WriteAndTrackFasta, gather_at_rank, get_ident)


GATHER_MIN_MATCHES=3
F_IDENT_THRESHOLD=0.1
F_MAJOR_THRESHOLD=0.2


from enum import Enum
class ContigInfo(Enum):
    "enum for contig information"
    CLEAN = 1                             # proven clean
    DIRTY = 2                             # proven dirty
    NO_IDENT = 3                          # no identified hashes in contig
    NO_HASH = 4                           # no hashes in this contig


def kb(bp):
    return int(bp/1000)


class ContigReport(object):
    # output in addition: genomefile, genome taxonomy, match_rank
    def __init__(self, contig_name, bp, info, reason, lineage, total_hashes,
                 ident_hashes, match_hashes):
        self.contig_name = contig_name
        self.bp = bp
        self.info = info
        self.reason = reason
        self.lineage = lineage
        self.total_hashes = total_hashes
        self.ident_hashes = ident_hashes
        self.match_hashes = match_hashes


def do_gather_breakdown(minhash, lca_db, lin_db, min_matches, genome_lineage,
                        match_rank, report_fp,
                        reverse_flag=False, display_fp=False):
    "Report all gather matches to report_fp; return first match count."
    threshold_percent = min_matches / len(minhash)

    # do the gather:
    first_count = None
    first_match_under_fp = False
    for lin, count in gather_at_rank(minhash, lca_db, lin_db, match_rank):
        if count <= min_matches:
            if display_fp and first_match_under_fp:
                first_match_under_fp = True
                print('  --------- (likely false positives below this line) ---------', file=report_fp)
            elif not display_fp:
                break

        if not first_count:               # set first_count once
            first_count = count

        assert lin[-1].rank == match_rank

        # display matches --
        f_count = count / len(minhash)
        lineage_display = lin[-1].name

        matches_lineage = True
        if not utils.is_lineage_match(genome_lineage, lin, match_rank):
            matches_lineage = False

        # should we flag this match?
        if not matches_lineage and not reverse_flag:
            lineage_display = '!! ' + lineage_display
        elif matches_lineage and reverse_flag:
            lineage_display = '!! ' + lineage_display

        print(f'  {f_count*100:.2f}% - to {lineage_display}', file=report_fp)

    if first_match_under_fp:
        print(f'** note: matches under {threshold_percent*100:.3f}% may be false positives', file=report_fp)

    return first_count


def create_empty_output(genome, comment, summary, report, contig_report,
                        clean, dirty,
                        f_major="", f_ident="",
                        provided_lin="", guessed_lineage=""):
    "Create empty output from early exit, so snakemake doesn't complain."
    if summary:
        with open(summary, 'wt') as fp:
            w = csv.writer(fp)
            if guessed_lineage:
                guessed_lineage = sourmash.lca.display_lineage(guessed_lineage)
            if provided_lin:
                provided_lin = sourmash.lca.display_lineage(provided_lin)

            row = [genome] + ["", f_major, f_ident] + [""]*14 + \
               [guessed_lineage, provided_lin, comment]
            w.writerow(row)
    if report:
        with open(report, 'wt') as fp:
            fp.write(comment)
    if contig_report:
        with open(contig_report, 'wt') as fp:
            pass
    open(clean, 'wt').close()
    open(dirty, 'wt').close()


def guess_tax_by_gather(entire_mh, lca_db, lin_db, match_rank, report_fp):
    "Guess likely taxonomy using gather."
    sum_ident = 0
    first_lin = ()
    first_count = 0

    for lin, count in gather_at_rank(entire_mh, lca_db, lin_db, match_rank):
        if count >= GATHER_MIN_MATCHES:
            # record the first lineage we come across as likely lineage.
            if not first_lin:
                first_lin = lin
                first_count = count

        sum_ident += count

    f_ident = sum_ident / len(entire_mh)
    f_major = first_count / sum_ident

    return first_lin, f_ident, f_major


class ContigsDecontaminator(object):
    GATHER_THRESHOLD = GATHER_MIN_MATCHES

    def __init__(self, genome_lineage, match_rank, empty_mh, lca_db, lin_db):
        self.genome_lineage = genome_lineage
        self.empty_mh = empty_mh
        self.match_rank = match_rank
        self.lca_db = lca_db
        self.lin_db = lin_db

        self.n_reason_1 = 0
        self.n_reason_2 = 0
        self.n_reason_3 = 0
        self.missed_n = 0
        self.missed_bp = 0
        self.noident_n = 0
        self.noident_bp = 0

        self.clean_out = None
        self.dirty_out = None

        self.contig_reports = {}

    def set_clean_filename(self, filename):
        clean_fp = gzip.open(filename, 'wt')
        self.clean_out = WriteAndTrackFasta(clean_fp, self.empty_mh)

    def set_dirty_filename(self, filename):
        dirty_fp = gzip.open(filename, 'wt')
        self.dirty_out = WriteAndTrackFasta(dirty_fp, self.empty_mh)

    def clean_contigs(self, screed_iter, report_fp, no_write=False):
        for n, record in enumerate(screed_iter):
            # make a new minhash and start examining it.
            mh = self.empty_mh.copy_and_clear()
            mh.add_sequence(record.sequence, force=True)

            clean_flag = ContigInfo.NO_HASH
            ctg_lin = ""
            hash_cnt = 0
            reason = 0

            if mh and len(mh) >= self.GATHER_THRESHOLD:
                clean_flag, ctg_lin, hash_cnt = self.check_gather(record, mh, report_fp)
                if clean_flag == ContigInfo.DIRTY: reason = 1

            # track things
            if clean_flag == ContigInfo.NO_HASH:
                self.missed_n += 1
                self.missed_bp += len(record.sequence)
            elif clean_flag == ContigInfo.NO_IDENT:
                self.noident_n += 1
                self.noident_bp += len(record.sequence)

            # write out contigs -> clean or dirty files.
            if clean_flag != ContigInfo.DIRTY:   # non-dirty => clean
                if self.clean_out:
                    self.clean_out.write(record, no_write=no_write)
            else:
                assert clean_flag == ContigInfo.DIRTY
                if self.dirty_out:
                    self.dirty_out.write(record, no_write=no_write)

            hash_ident_cnt = 0
            for hashval in mh.get_mins():
                if hashval in self.lca_db.hashval_to_idx:
                    hash_ident_cnt += 1
            ctg_rep = ContigReport(record.name, len(record.sequence),
                                   clean_flag, reason, ctg_lin,
                                   len(mh), hash_ident_cnt, hash_cnt)
            assert record.name not in self.contig_reports
            self.contig_reports[record.name] = ctg_rep

        # END contig loop

    def check_gather(self, record, contig_mh, report_fp):
        """
        Does this contig have a gather match that is outside the given rank?

        This method discovers chunks of sequence that belong to other
        lineages.
        """
        if not contig_mh:
            return ContigInfo.NO_HASH, "", 0

        good_count = 0
        for lin, count in gather_at_rank(contig_mh, self.lca_db,
                                         self.lin_db, self.match_rank):
            if utils.is_lineage_match(lin, self.genome_lineage,
                                      self.match_rank):
                good_count += count

        if good_count >= self.GATHER_THRESHOLD:
            return ContigInfo.CLEAN, "", good_count

        threshold_bp = contig_mh.scaled * self.GATHER_THRESHOLD
        results = self.lca_db.gather(sourmash.SourmashSignature(contig_mh),
                                     threshold_bp=threshold_bp)

        if not results:
            return ContigInfo.NO_IDENT, "", 0

        score = results[0][0]
        match = results[0][1]

        if score < F_IDENT_THRESHOLD:
            return ContigInfo.NO_IDENT, "", 0

        # get identity
        match_ident = get_ident(match)
        # get lineage
        contig_lineage = self.lin_db.ident_to_lineage[match_ident]
        common_hashcount = contig_mh.count_common(match.minhash)
        common_kb = common_hashcount * contig_mh.scaled / 1000
        # if it matched outside rank, => dirty.
        if utils.is_lineage_match(self.genome_lineage, contig_lineage,
                                  self.match_rank):
            clean = ContigInfo.CLEAN
        else:
            clean = ContigInfo.DIRTY
            self.n_reason_1 += 1

            print(f'---- contig {record.name} ({len(record.sequence)/1000:.0f} kb)\n', file=report_fp)
            print(f'contig dirty, REASON 1 - gather matches to lineage outside of genome\'s {self.match_rank}', file=report_fp)
            print(f'   contig gather yields match of {common_kb:.0f} kb to {pretty_print_lineage2(contig_lineage, self.match_rank)}', file=report_fp)
            print(f'   vs genome lineage of {pretty_print_lineage2(self.genome_lineage, self.match_rank)}', file=report_fp)
            print('', file=report_fp)

        return clean, contig_lineage, common_hashcount

    def _report_lca_summary(self, report_fp, ctg_tax_assign, ctg_assign,
                            f_match, f_ident):
        scaled = self.empty_mh.scaled

        ctg_counts = Counter()
        for hashval, lineages in ctg_assign.items():
            for lineage in lineages:
                ctg_counts[lineage] += 1

        print(f'\nf_ident {f_ident:.2f} / f_match {f_match:.2f}', file=report_fp)
        print(f'\n** hashval lca counts', file=report_fp)
        for lin, count in ctg_tax_assign.most_common():
            print(f'   {count*scaled/1000:.0f} kb {pretty_print_lineage2(lin, self.match_rank)}', file=report_fp)
        print(f'\n** hashval lineage counts - {len(ctg_assign)}', file=report_fp)
        for lin, count in ctg_counts.most_common():
            print(f'   {count*scaled/1000:.0f} kb {pretty_print_lineage2(lin, self.match_rank)}', file=report_fp)
        print('', file=report_fp)

# END CLASS


def choose_genome_lineage(guessed_genome_lineage, provided_lineage, match_rank,
                          f_ident, f_major, report):

    comment = ""
    genome_lineage = None

    if provided_lineage:
        if utils.is_lineage_match(provided_lineage, guessed_genome_lineage, match_rank):
            report(f'(provided lineage agrees with k-mer classification at {match_rank} level)')
        elif guessed_genome_lineage:
            report(f'(provided lineage disagrees with k-mer classification at or above {match_rank} level)')
        else:
            pass

        genome_lineage = utils.pop_to_rank(provided_lineage, match_rank)
        report(f'\nUsing provided lineage as genome lineage.')
    else:
        if f_ident < F_IDENT_THRESHOLD:
            report(f'** ERROR: fraction of total identified hashes (f_ident) < {F_IDENT_THRESHOLD*100:.0f}%.')
            comment = f"too few identifiable hashes; f_ident < {F_IDENT_THRESHOLD*100:.0f}%. provide a lineage for this genome."
        elif f_major < F_MAJOR_THRESHOLD:
            report(f'** ERROR: fraction of identified hashes in major lineage (f_major) < {F_MAJOR_THRESHOLD*100:.0f}%.')
            comment = f"too few hashes in major lineage; f_major < {F_MAJOR_THRESHOLD*100:.0f}%. provide a lineage for this genome."
        else:
            genome_lineage = utils.pop_to_rank(guessed_genome_lineage, match_rank)
            report(f'Using majority gather lineage as genome lineage.')

    return genome_lineage, comment


###

def main(args):
    "Main entry point for scripting. Use cmdline for command line entry."
    genomebase = os.path.basename(args.genome)
    match_rank = args.match_rank

    tax_assign, _ = load_taxonomy_assignments(args.lineages_csv,
                                              start_column=3)
    print(f'loaded {len(tax_assign)} tax assignments.')

    with open(args.matches_sig, 'rt') as fp:
        siglist = list(sourmash.load_signatures(fp))

    if not siglist:
        print('no matches for this genome, exiting.')
        comment = "no matches to this genome were found in the database; nothing to do"
        create_empty_output(genomebase, comment, args.summary,
                            args.report, args.contig_report,
                            args.clean, args.dirty)
        return 0

    report_fp = open(args.report, 'wt')
    def report(*args):
        print(*args)
        print(*args, file=report_fp)

    # construct a template minhash object that we can use to create new 'uns
    empty_mh = siglist[0].minhash.copy_and_clear()
    ksize = empty_mh.ksize
    scaled = empty_mh.scaled
    moltype = empty_mh.moltype

    # create empty LCA database to populate...
    lca_db = LCA_Database(ksize=ksize, scaled=scaled, moltype=moltype)
    lin_db = LineageDB()

    # ...with specific matches.
    for ss in siglist:
        ident = get_ident(ss)
        lineage = tax_assign[ident]

        lca_db.insert(ss, ident=ident)
        lin_db.insert(ident, lineage)

    print(f'loaded {len(siglist)} signatures & created LCA Database')

    report(f'charcoal version: v{version}')
    report(f'match_rank: {match_rank} / scaled: {scaled} / ksize: {ksize}')
    report('')
    report(f'genome: {genomebase}')

    print(f'pass 1: reading contigs from {genomebase}')
    entire_mh = empty_mh.copy_and_clear()
    total_bp = 0
    for n, record in enumerate(screed.open(args.genome)):
        entire_mh.add_sequence(record.sequence, force=True)
        total_bp += len(record.sequence)
    n_contigs = n + 1

    report(f'genome has {n_contigs} contigs in {total_bp/1000:.1f}kb')
    report(f'{len(entire_mh)} hashes total.')
    report('')

    # calculate lineage from majority vote on LCA
    guessed_genome_lineage, f_major, f_ident = \
         guess_tax_by_gather(entire_mh, lca_db, lin_db, match_rank, report_fp)

    report(f'Gather classification on this genome yields: {pretty_print_lineage(guessed_genome_lineage)}')

    if f_major == 1.0 and f_ident == 1.0:
        comment = "All genome hashes belong to one lineage! Nothing to do."
        report(comment)
        create_empty_output(genomebase, comment, args.summary,
                            None, args.contig_report,
                            args.clean, args.dirty,
                            guessed_lineage=guessed_genome_lineage,
                            f_ident=f_ident, f_major=f_major)
        return 0

    # did we get a passed-in lineage assignment?
    provided_lin = ""
    if args.lineage and args.lineage != 'NA':
        provided_lin = args.lineage.split(';')
        provided_lin = [ LineagePair(rank, name) for (rank, name) in zip(sourmash.lca.taxlist(), provided_lin) if name.strip() ]
        report(f'Provided lineage from command line:\n   {sourmash.lca.display_lineage(provided_lin)}')

    # choose between the lineages
    genome_lineage, comment = choose_genome_lineage(guessed_genome_lineage,
                                                    provided_lin,
                                                    match_rank,
                                                    f_ident, f_major,
                                                    report)

    if comment: # failure to get a good lineage assignment? exit early.
        report(f'** Please provide a lineage for this genome.')

        # ...unless we force.
        if args.force:
            print('--force requested, so continuing despite this.')
        else:
            create_empty_output(genomebase, comment, args.summary,
                                None, args.contig_report,
                                args.clean, args.dirty,
                                provided_lin=provided_lin,
                                guessed_lineage=guessed_genome_lineage,
                                f_ident=f_ident, f_major=f_major)
            return 0


    # is match_rank lower than genome lineage? if so, raise it.
    lineage_ranks = [ x.rank for x in genome_lineage ]
    if match_rank not in lineage_ranks:
        report(f'** NOTE: lineage rank is {lineage_ranks[-1]}; pulling match rank back to that.')
        match_rank = lineage_ranks[-1]

    report(f'\nFull lineage being used for contamination analysis:')
    report(f'   {sourmash.lca.display_lineage(genome_lineage)}')

    cleaner = ContigsDecontaminator(genome_lineage, match_rank,
                                    empty_mh, lca_db, lin_db)
    cleaner.set_clean_filename(args.clean)
    cleaner.set_dirty_filename(args.dirty)

    print('')
    print(f'pass 2: reading contigs from {genomebase}')
    print(f'\n**\n** walking through contigs:\n**\n', file=report_fp)

    # do the cleaning
    screed_iter = screed.open(args.genome)
    cleaner.clean_contigs(screed_iter, report_fp)
    dirty_mh = cleaner.dirty_out.minhash

    fail = False
    for lin, count in gather_at_rank(dirty_mh, lca_db, lin_db, match_rank):
        if count >= GATHER_MIN_MATCHES and utils.is_lineage_match(genome_lineage,
                                                           lin, match_rank):
            fail = True

    if fail and 0:
        print(f'\nbreakdown of dirty contigs w/gather:', file=report_fp)
        do_gather_breakdown(dirty_mh, lca_db, lin_db,
                            GATHER_MIN_MATCHES,
                            genome_lineage, match_rank,
                            report_fp)

        comment = "cannot remove contamination without removing clean sequence; punting"
        report(comment)
        create_empty_output(genomebase, comment, args.summary,
                            args.report, args.contig_report, args.clean,
                            args.dirty)
        sys.exit(0)

    # recover information from cleaner object
    clean_n = cleaner.clean_out.n
    clean_bp = cleaner.clean_out.bp
    dirty_n = cleaner.dirty_out.n
    dirty_bp = cleaner.dirty_out.bp
    missed_n = cleaner.missed_n
    missed_bp = cleaner.missed_bp
    noident_n = cleaner.noident_n
    noident_bp = cleaner.noident_bp

    n_reason_1 = cleaner.n_reason_1
    n_reason_2 = cleaner.n_reason_2
    n_reason_3 = cleaner.n_reason_3

    assert n_reason_1 + n_reason_2 + n_reason_3 == dirty_n

    # do some reporting.
    report('--------------')

    report('')
    report(f'Summary report:')
    report('')
    report(f'genome {match_rank} is {genome_lineage[-1].name}.')
    report('')

    report(f'kept {clean_n} contigs containing {int(clean_bp/1000)} kb.')
    report(f'removed {dirty_n} contigs containing {int(dirty_bp/1000)} kb.')
    report(f'{missed_n} contigs ({int(missed_bp/1000)} kb total) had no hashes, and were counted as clean')

    # look at what our database says about remaining contamination,
    # across all "clean" contigs.

    # report gather breakdown of dirty signature
    print(f'\nbreakdown of dirty sequence w/gather ({dirty_n} contigs, {kb(dirty_bp)} kb):', file=report_fp)

    dirty_mh = cleaner.dirty_out.minhash
    first_match = None
    if dirty_mh:
        first_match = do_gather_breakdown(dirty_mh, lca_db, lin_db,
                            GATHER_MIN_MATCHES,
                            genome_lineage, match_rank,
                            report_fp, reverse_flag=True)
    if not first_match:
        print(' ** no matches to dirty sequence **', file=report_fp)

    # report gather breakdown of clean signature
    print(f'\nbreakdown of clean sequence w/gather ({clean_n} contigs, {kb(clean_bp)} kb):', file=report_fp)

    clean_mh = cleaner.clean_out.minhash
    first_match = None
    if clean_mh:
        first_match = do_gather_breakdown(clean_mh, lca_db, lin_db,
                                          GATHER_MIN_MATCHES,
                                          genome_lineage, match_rank,
                                          report_fp)

    if not first_match:
        print(' ** no matches to clean sequence **', file=report_fp)

    report('')
    report('Please have a nice day!')

    ### output:

    # get genome size and match lineage of primary clean match
    nearest_size = 0
    match_lineage = ""
    ratio = 0.0
    if first_match:
        nearest_size = first_match * scaled
        ident = '**'
        match_lineage = ""
        ratio = 0

    # write out a one line summary:
    if args.summary:
        comment = ""

        with open(args.summary, 'wt') as fp:
            full_lineage = sourmash.lca.display_lineage(match_lineage)
            short_lineage = pretty_print_lineage(genome_lineage)
            f_removed = dirty_bp / (dirty_bp + clean_bp)
            w = csv.writer(fp)
            w.writerow([genomebase, short_lineage,
                        f_major, f_ident, f_removed,
                        n_reason_1, n_reason_2, n_reason_3,
                        nearest_size, ratio, clean_bp,
                        clean_n, dirty_n, dirty_bp,
                        missed_n, missed_bp,
                        noident_n, noident_bp,
                        full_lineage,
                        sourmash.lca.display_lineage(provided_lin),
                        comment])

    if args.contig_report:
        with open(args.contig_report, 'wt') as fp:
            w = csv.writer(fp)
            w.writerow(['genomefile',
                           'genometax',
                           'match_rank',
                           'contig_name',
                           'bp',
                           'decision',
                           'reason',
                           'contigtax',
                           'total_hashes',
                           'ident_hashes',
                           'match_hashes'])
            for rep in cleaner.contig_reports.values():
                w.writerow([args.genome,
                            sourmash.lca.display_lineage(genome_lineage),
                            match_rank,
                            rep.contig_name,
                            rep.bp,
                            rep.info,
                            rep.reason,
                            sourmash.lca.display_lineage(rep.lineage),
                            rep.total_hashes,
                            rep.ident_hashes,
                            rep.match_hashes])

    return 0


def cmdline(sys_args):
    "Command line entry point w/argparse action."
    p = argparse.ArgumentParser(sys_args)
    p.add_argument('--genome', help='genome file', required=True)
    p.add_argument('--lineages_csv', help='lineage spreadsheet', required=True)
    p.add_argument('--matches_sig', help='all relevant matches', required=True)
    p.add_argument('--clean', help='cleaned contigs', required=True)
    p.add_argument('--dirty', help='dirty contigs', required=True)
    p.add_argument('--report', help='report output', required=True)
    p.add_argument('--summary', help='CSV one line output')
    p.add_argument('--force', help='continue past survivable errors',
                   action='store_true')

    p.add_argument('--lineage', help=';-separated lineage down to genus level',
                   default='NA')        # default is str NA
    p.add_argument('--match-rank', help='rank below which matches are _not_ contaminants', default='genus')
    p.add_argument('--contig-report', help='contig report (CSV)')
    args = p.parse_args()

    main(args)


# execute this, when run with `python -m`.
if __name__ == '__main__':
    returncode = cmdline(sys.argv[1:])
    sys.exit(0)
