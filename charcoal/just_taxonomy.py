#! /usr/bin/env python
"""
Remove bad contigs based solely on taxonomy.

CTB TODO:
* optionally eliminate contigs with no taxonomy
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


GATHER_MIN_MATCHES=3


def get_idents_for_hashval(lca_db, hashval):
    "Get the identifiers associated with this hashval."
    idx_list = lca_db.hashval_to_idx.get(hashval, [])
    for idx in idx_list:
        ident = lca_db.idx_to_ident[idx]
        yield ident


def gather_assignments(hashvals, rank, dblist, ldb):
    """
    Collect lineage assignments from across all the databases for all the
    hashvals.
    """
    assignments = defaultdict(set)
    for hashval in hashvals:
        for lca_db in dblist:
            lineages = set()
            for ident in get_idents_for_hashval(lca_db, hashval):
                lineage = ldb.ident_to_lineage[ident]

                if rank:
                    lineage = utils.pop_to_rank(lineage, rank)
                assignments[hashval].add(lineage)

    return assignments


def count_lca_for_assignments(assignments):
    """
    For each hashval, count the LCA across its assignments.
    """
    counts = Counter()
    for hashval in assignments:

        # for each list of tuple_info [(rank, name), ...] build
        # a tree that lets us discover lowest-common-ancestor.
        lineages = assignments[hashval]
        tree = sourmash.lca.build_tree(lineages)

        # now find either a leaf or the first node with multiple
        # children; that's our lowest-common-ancestor node.
        lca, reason = sourmash.lca.find_lca(tree)
        counts[lca] += 1

    return counts


def pretty_print_lineage(lin):
    "Nice output names for lineages."
    if not lin:
        return f'** no assignment **'
    elif lin[-1].rank == 'strain':
        strain = lin[-1].name
        return f'{strain}'
    elif lin[-1].rank == 'species':
        species = lin[-1].name
        return f'{species}'
    else:
        return f'{lin[-1].rank} {lin[-1].name}'


def get_ident(sig):
    "Hack and slash identifiers."
    ident = sig.name()
    ident = ident.split()[0]
    ident = ident.split('.')[0]
    return ident


def check_gather(record, contig_mh, genome_lineage, match_rank,
                 lca_db, lineage_db, report_fp):
    "Does this contig have a gather match that is outside the given rank?"
    threshold_bp = contig_mh.scaled*GATHER_MIN_MATCHES
    results = lca_db.gather(sourmash.SourmashSignature(contig_mh),
                            threshold_bp=threshold_bp)

    if not results:
        return True

    match = results[0][1]

    # get identity
    match_ident = get_ident(match)
    # get lineage
    contig_lineage = lineage_db.ident_to_lineage[match_ident]

    # if it matched outside rank, => dirty.
    clean = True
    if not utils.is_lineage_match(genome_lineage, contig_lineage,
                                  match_rank):
        clean=False
        common_kb = contig_mh.count_common(match.minhash) * contig_mh.scaled / 1000

        print(f'---- contig {record.name} ({len(record.sequence)/1000:.0f} kb)', file=report_fp)
        print(f'contig dirty, REASON 3 - gather matches to lineage outside of genome\'s {match_rank}\n   gather yields match of {common_kb:.0f} kb to {pretty_print_lineage(contig_lineage)}',
              file=report_fp)
        print('', file=report_fp)

    return clean


def check_lca(record, contig_mh, genome_lineage, match_rank,
              lca_db, lin_db, report_fp):
    "Does this contig have any hashes with LCA outside the given rank?"
    clean = True
    reason = 0

    # get _all_ of the hash taxonomy assignments for this contig
    ctg_assign = gather_assignments(contig_mh.get_mins(), None,
                                    [lca_db], lin_db)

    ctg_tax_assign = count_lca_for_assignments(ctg_assign)
    if not ctg_tax_assign:
        return True, ""

    # get top assignment for contig.
    ctg_lin, lin_count = next(iter(ctg_tax_assign.most_common()))

    # assignment outside of genus? dirty!
    ok_ranks = set()
    passed_rank = False
    for rank in sourmash.lca.taxlist():
        if rank == match_rank:
            passed_rank = True

        if passed_rank:
            ok_ranks.add(rank)

    if (not ctg_lin) or (ctg_lin[-1].rank not in ok_ranks):
        bad_rank = "(root)"
        if ctg_lin:
            bad_rank = ctg_lin[-1].rank
        clean = False
        reason = 1
        print(f'\n---- contig {record.name} ({len(record.sequence)/1000:.0f} kb)', file=report_fp)
        print(f'contig dirty, REASON 1 - contig LCA is above {match_rank}\nlca rank is {bad_rank}',
              file=report_fp)
    elif not utils.is_lineage_match(genome_lineage, ctg_lin, match_rank):
        clean = False
        reason = 2
        print('', file=report_fp)
        print(f'---- contig {record.name} ({len(record.sequence)/1000:.0f} kb)', file=report_fp)
        print(f'contig dirty, REASON 2 - contig lineage is not a match to genome\'s {match_rank}\nlineage is {pretty_print_lineage(ctg_lin)}',
              file=report_fp)

    # summary reporting --
    if not clean:
        report_lca_summary(report_fp, ctg_tax_assign, ctg_assign,
                           contig_mh.scaled)

    return clean, reason


def report_lca_summary(report_fp, ctg_tax_assign, ctg_assign, scaled):
    ctg_counts = Counter()
    for hashval, lineages in ctg_assign.items():
        for lineage in lineages:
            ctg_counts[lineage] += 1

    print(f'\n** hashval lca counts', file=report_fp)
    for lin, count in ctg_tax_assign.most_common():
        print(f'   {count*scaled/1000:.0f} kb {pretty_print_lineage(lin)}', file=report_fp)
    print(f'\n** hashval lineage counts - {len(ctg_assign)}', file=report_fp)
    for lin, count in ctg_counts.most_common():
        print(f'   {count*scaled/1000:.0f} kb {pretty_print_lineage(lin)}', file=report_fp)
    print('', file=report_fp)


class WriteAndTrackFasta(object):
    def __init__(self, outfp, mh_ex):
        self.minhash = mh_ex.copy_and_clear()
        self.outfp = outfp
        self.n = 0
        self.bp = 0

    def write(self, record):
        self.outfp.write(f'>{record.name}\n{record.sequence}\n')
        self.minhash.add_sequence(record.sequence, force=True)
        self.n += 1
        self.bp += len(record.sequence)

    def close(self):
        self.outfp.close()


def do_gather_breakdown(minhash, lca_db, report_fp):
    "Report all gather matches to report_fp; return first match sig."
    import copy
    minhash = copy.copy(minhash)
    query_sig = sourmash.SourmashSignature(minhash)

    threshold_percent = GATHER_MIN_MATCHES  / len(minhash)

    # do the gather:
    first_match = None
    first_match_under_fp = False
    while 1:
        results = lca_db.gather(query_sig, threshold_bp=0)
        if not results:
            break

        (match, match_sig, _) = results[0]
        if not first_match:               # set first_match once
            first_match = match_sig

        if match <= threshold_percent and not first_match_under_fp:
            first_match_under_fp = True
            print('  --------- (likely false positives below this line) ---------', file=report_fp)

        print(f'  {match*100:.2f}% - to {match_sig.name()}', file=report_fp)
        minhash.remove_many(match_sig.minhash.get_mins())
        query_sig = sourmash.SourmashSignature(minhash)

    if first_match_under_fp:
        print(f'** note: matches under {threshold_percent*100:.3f}% may be false positives', file=report_fp)

    return first_match


def create_empty_output(genome, comment, summary, report, clean, dirty,
                        f_major="", f_ident="",
                        provided_lin="", lca_lineage=""):
    if summary:
        with open(summary, 'wt') as fp:
            w = csv.writer(fp)
            if lca_lineage:
                lca_lineage = sourmash.lca.display_lineage(lca_lineage)
            if provided_lin:
                provided_lin = sourmash.lca.display_lineage(provided_lin)

            row = [genome] + ["", f_major, f_ident] + [""]*12 + \
               [lca_lineage, provided_lin, comment]
            w.writerow(row)
    if report:
        with open(report, 'wt') as fp:
            fp.write(comment)
    open(clean, 'wt').close()
    open(dirty, 'wt').close()


def get_majority_lca_at_rank(entire_mh, lca_db, lin_db, rank, report_fp):
    # get all of the hash taxonomy assignments for this genome
    hash_assign = gather_assignments(entire_mh.get_mins(), rank,
                                     [lca_db], lin_db)

    # count them and find major
    counts = Counter()
    identified_counts = 0
    for hashval, lineages in hash_assign.items():
        if lineages:
            identified_counts += 1
            for lineage in lineages:
                lineage = utils.pop_to_rank(lineage, 'genus')
                counts[lineage] += 1

    genome_lineage, count = next(iter(counts.most_common()))

    f_ident = identified_counts / len(entire_mh)
    f_major = count / identified_counts
    total_counts = len(hash_assign)

    # report everything...

    print(f'{f_ident*100:.1f}% of total hashes identified.', file=report_fp)
    print(f'{f_major*100:.1f}% of identified hashes match to {pretty_print_lineage(genome_lineage)}', file=report_fp)
    print(f'({identified_counts} identified hashes, {count} in most common)', file=report_fp)
    if f_major < 0.8:
        print(f'** WARNING ** majority lineage is less than 80% of assigned lineages. Beware!', file=report_fp)

    print(f'\n** hashval lineage counts for genome - {total_counts} => {total_counts*entire_mh.scaled/1000:.0f} kb', file=report_fp)
    for lin, count in counts.most_common():
        print(f'   {count*entire_mh.scaled/1000:.0f} kb {pretty_print_lineage(lin)}', file=report_fp)
    print('', file=report_fp)

    return genome_lineage, f_major, f_ident


def main(args):
    p = argparse.ArgumentParser(args)
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
    args = p.parse_args()

    genomebase = os.path.basename(args.genome)

    tax_assign, _ = load_taxonomy_assignments(args.lineages_csv,
                                              start_column=3)
    print(f'loaded {len(tax_assign)} tax assignments.')

    with open(args.matches_sig, 'rt') as fp:
        siglist = list(sourmash.load_signatures(fp))

    if not siglist:
        print('no matches for this genome, exiting.')
        comment = "no matches to this genome were found in the database; nothing to do"
        create_empty_output(genomebase, comment, args.summary,
                            args.report, args.clean, args.dirty)
        return 0

    report_fp = open(args.report, 'wt')
    def report(*args):
        print(*args)
        print(*args, file=report_fp)

    # construct a template minhash object that we can use to create new 'uns
    empty_mh = siglist[0].minhash.copy_and_clear()
    ksize = empty_mh.ksize
    scaled = empty_mh.scaled

    # create empty LCA database to populate...
    lca_db = LCA_Database(ksize=ksize, scaled=scaled)
    lin_db = LineageDB()

    # ...with specific matches.
    for ss in siglist:
        ident = get_ident(ss)
        lineage = tax_assign[ident]

        lca_db.insert(ss, ident=ident)
        lin_db.insert(ident, lineage)

    print(f'loaded {len(siglist)} signatures & created LCA Database')

    print(f'pass 1: reading contigs from {genomebase}')
    entire_mh = empty_mh.copy_and_clear()
    for n, record in enumerate(screed.open(args.genome)):
        entire_mh.add_sequence(record.sequence, force=True)

    # calculate lineage from majority vote on LCA
    lca_genome_lineage, f_major, f_ident = \
         get_majority_lca_at_rank(entire_mh, lca_db, lin_db, 'genus',
                                  report_fp)

    report(f'K-mer classification on this genome yields: {pretty_print_lineage(lca_genome_lineage)}')

    # did we get a passed-in lineage assignment?
    provided_lin = ""
    if args.lineage and args.lineage != 'NA':
        provided_lin = args.lineage.split(';')
        provided_lin = [ LineagePair(rank, name) for (rank, name) in zip(sourmash.lca.taxlist(), provided_lin) ]
        report(f'Provided lineage from command line:\n   {sourmash.lca.display_lineage(provided_lin)}')

        if utils.is_lineage_match(provided_lin, lca_genome_lineage, 'genus'):
            report(f'(provided lineage agrees with k-mer classification at genus level)')
        else:
            report(f'(provided lineage disagrees with k-mer classification at or above genus level)')

        genome_lineage = utils.pop_to_rank(provided_lin, 'genus')
        report(f'\nUsing provided lineage as genome lineage.')
    else:
        fail = False
        if f_ident < 0.1:
            print(f'** ERROR: fraction of total identified hashes (f_ident) < 10%.')
            print(f'** Please provide a lineage for this genome.')
            print(f'** ERROR: fraction of identified hashes f_major < 20%.',
                  file=report_fp)
            print(f'** Please provide a lineage for this genome.',
                  file=report_fp)
            comment = "too few identifiable hashes; < 10%. provide a lineage for this genome."
            fail = True
        elif f_major < 0.2:
            print(f'** ERROR: fraction of identified hashes in major lineage (f_major) < 20%.')
            print(f'** Please provide a lineage for this genome.')
            print(f'** ERROR: fraction of identified hashes f_major < 20%.',
                  file=report_fp)
            print(f'** Please provide a lineage for this genome.',
                  file=report_fp)
            comment = "too few hashes in major lineage; < 20%. provide a lineage for this genome."
            fail = True

        if fail:
            if args.force:
                print('--force requested, so continuing despite this.')
            else:
                create_empty_output(genomebase, comment, args.summary,
                                    None, args.clean, args.dirty,
                                    provided_lin=provided_lin,
                                    lca_lineage=lca_genome_lineage,
                                    f_ident=f_ident, f_major=f_major)
                
                return 0

        genome_lineage = lca_genome_lineage
        report(f'Using LCA majority lineage as genome lineage.')


    # is match_rank lower than genome lineage? if so, raise it.
    match_rank = args.match_rank
    lineage_ranks = [ x.rank for x in genome_lineage ]
    if match_rank not in lineage_ranks:
        report(f'** NOTE: lineage rank is {lineage_ranks[-1]}; pulling match rank back to that.')
        match_rank = lineage_ranks[-1]

    report(f'\nFull lineage being used for contamination analysis:')
    report(f'   {sourmash.lca.display_lineage(genome_lineage)}')

    # the output files are coming!
    clean_fp = gzip.open(args.clean, 'wt')
    clean_out = WriteAndTrackFasta(clean_fp, empty_mh)
    dirty_fp = gzip.open(args.dirty, 'wt')
    dirty_out = WriteAndTrackFasta(dirty_fp, empty_mh)

    missed_n = 0
    missed_bp = 0

    # now, find bad contigs.
    n_reason_1 = 0
    n_reason_2 = 0
    n_reason_3 = 0

    print('')
    print(f'pass 2: reading contigs from {genomebase}')
    print(f'\n**\n** walking through contigs:\n**\n', file=report_fp)
    for n, record in enumerate(screed.open(args.genome)):
        # make a new minhash and start examining it.
        mh = empty_mh.copy_and_clear()
        mh.add_sequence(record.sequence, force=True)

        clean = True               # default to clean
        if not mh:                 # no hashes?
            missed_n += 1
            missed_bp += len(record.sequence)

        if mh and len(mh) >= GATHER_MIN_MATCHES: # CTB: don't hard code.
            clean = check_gather(record, mh, genome_lineage, match_rank,
                                 lca_db, lin_db, report_fp)
            if not clean:
                n_reason_3 += 1

        # did we find a dirty contig in step 1? if NOT, go into LCA style
        # approaches.
        if mh and clean:
            clean, reason = check_lca(record, mh, genome_lineage, match_rank,
                                      lca_db, lin_db, report_fp)
            if not clean:
                if reason == 1:
                    n_reason_1 += 1
                elif reason == 2:
                    n_reason_2 += 1
                else:
                    assert 0, "unknown dirty reason code"

        # write out contigs -> clean or dirty files.
        if clean:
            clean_out.write(record)
        else:
            dirty_out.write(record)

    # END contig loop

    clean_n = clean_out.n
    clean_bp = clean_out.bp
    dirty_n = dirty_out.n
    dirty_bp = dirty_out.bp

    assert n_reason_1 + n_reason_2 + n_reason_3 == dirty_n

    # do some reporting.
    print('--------------', file=report_fp)
    report(f'kept {clean_n} contigs containing {int(clean_bp/1000)} kb.')
    report(f'removed {dirty_n} contigs containing {int(dirty_bp/1000)} kb.')
    report(f'{missed_n} contigs ({int(missed_bp/1000)} kb total) had no hashes, and counted as clean')

    # look at what our database says about remaining contamination,
    # across all "clean" contigs. (CTB: Need to dig into this more to figure
    # out exactly why we still have any :)
    # CTB: add breakdown of dirty contigs?

    # report gather breakdown of clean signature
    print(f'\nbreakdown of clean contigs w/gather:', file=report_fp)

    clean_mh = clean_out.minhash
    first_match = None
    if clean_mh:
        first_match = do_gather_breakdown(clean_mh, lca_db, report_fp)

    if not first_match:
        print(' ** no matches **', file=report_fp)

    # get genome size and match lineage of primary match
    nearest_size = 0
    match_lineage = ""
    ratio = 0.0
    if first_match:
        nearest_size = len(first_match.minhash) * first_match.minhash.scaled
        ident = get_ident(first_match)
        match_lineage = lin_db.ident_to_lineage[ident]
        ratio = round(clean_bp / nearest_size, 2)

    # write out a one line summary?
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
                        full_lineage,
                        sourmash.lca.display_lineage(provided_lin),
                        comment])

    return 0


if __name__ == '__main__':
    returncode = main(sys.argv[1:])
    sys.exit(0)
