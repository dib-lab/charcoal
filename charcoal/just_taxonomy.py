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

import screed

import sourmash
from sourmash.lca.command_index import load_taxonomy_assignments
from sourmash.lca import LCA_Database, LineagePair

from . import utils
from . import lineage_db
from .lineage_db import LineageDB


def get_idents_for_hashval(lca_db, hashval):
    "Get the identifiers associated with this hashval."
    idx_list = lca_db.hashval_to_idx.get(hashval, [])
    for idx in idx_list:
        ident = lca_db.idx_to_ident[idx]
        yield ident


def gather_assignments(hashvals, rank, dblist, ldb):
    """
    Gather lineage assignments from across all the databases for all the
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
    if lin[-1].rank == 'strain':
        strain = lin[-1].name
        return f'{strain}'
    elif lin[-1].rank == 'species':
        species = lin[-1].name
        return f'{species}'
    elif not lin:
        return f'** no assignment **'
    else:
        return f'{lin[-1].rank} {lin[-1].name}'


def get_ident(sig):
    "Hack and slash identifiers."
    ident = sig.name()
    ident = ident.split()[0]
    ident = ident.split('.')[0]
    return ident


def check_gather(record, contig_mh, genome_lineage, lca_db, lineage_db, report_fp):
    threshold_bp = contig_mh.scaled*2
    results = lca_db.gather(sourmash.SourmashSignature(contig_mh))

    if not results:
        return True

    match = results[0][1]

    # get identitiy
    match_ident = get_ident(match)
    # get lineage
    contig_lineage = lineage_db.ident_to_lineage[match_ident]

    # if it matched outside genus, => dirty.
    clean = True
    if not utils.is_lineage_match(genome_lineage, contig_lineage, 'genus'):
        clean=False
        common_kb = contig_mh.count_common(match.minhash) * contig_mh.scaled / 1000

        print(f'---- contig {record.name} ({len(record.sequence)/1000:.0f} kb)', file=report_fp)
        print(f'contig dirty, REASON 3 - gather matches to lineage outside of genome\'s genus\n   gather yields match of {common_kb:.0f} kb to {pretty_print_lineage(contig_lineage)}',
              file=report_fp)
        print('', file=report_fp)

    return clean


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
    import copy
    minhash = copy.copy(minhash)
    query_sig = sourmash.SourmashSignature(minhash)

    # do the gather:
    first_match = None
    while 1:
        results = lca_db.gather(query_sig, threshold_bp=0)
        if not results:
            break

        (match, match_sig, _) = results[0]
        if not first_match:
            first_match = match_sig

        print(f'  {match*100:.3f}% - to {match_sig.name()}', file=report_fp)
        minhash.remove_many(match_sig.minhash.get_mins())
        query_sig = sourmash.SourmashSignature(minhash)

    return first_match


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--genome', help='genome file', required=True)
    p.add_argument('--lineages_csv', help='lineage spreadsheet', required=True)
    p.add_argument('--matches_sig', help='all relevant matches', required=True)
    p.add_argument('--clean', help='cleaned contigs', required=True)
    p.add_argument('--dirty', help='dirty contigs', required=True)
    p.add_argument('--report', help='report output', required=True)
    p.add_argument('--summary', help='CSV one line output')

    p.add_argument('--lineage', help=';-separated lineage down to genus level')
    args = p.parse_args()

    tax_assign, _ = load_taxonomy_assignments(args.lineages_csv,
                                              start_column=3)
    print(f'loaded {len(tax_assign)} tax assignments.')

    with open(args.matches_sig, 'rt') as fp:
        siglist = list(sourmash.load_signatures(fp))

    if not siglist:
        print('no matches for this genome, exiting.')
        comment = "no matches to this genome were found in the database"
        if args.summary:
            with open(args.summary, 'wt') as fp:
                w = csv.writer(fp)
                w.writerow([args.genome] + [""]*14 + [comment])
        open(args.report, 'wt').close()
        open(args.clean, 'wt').close()
        open(args.dirty, 'wt').close()
        sys.exit(0)

    # construct a template minhash object that we can use to create new 'uns
    empty_mh = siglist[0].minhash.copy_and_clear()
    ksize = empty_mh.ksize
    scaled = empty_mh.scaled

    # create empty LCA database to populate...
    lca_db = LCA_Database(ksize=ksize, scaled=scaled)
    ldb = LineageDB()

    # ...with specific matches.
    for ss in siglist:
        ident = get_ident(ss)
        lineage = tax_assign[ident] # could move functionality to LineageDB

        lca_db.insert(ss, ident=ident)
        ldb.insert(ident, lineage)

    print(f'loaded {len(siglist)} signatures & created LCA Database')

    print(f'pass 1: reading contigs from {args.genome}')
    entire_mh = empty_mh.copy_and_clear()
    for n, record in enumerate(screed.open(args.genome)):
        entire_mh.add_sequence(record.sequence, force=True)

    # get all of the hash taxonomy assignments for this genome
    hash_assign = gather_assignments(entire_mh.get_mins(), 'genus',
                                     [lca_db], ldb)

    # count them and find major
    counts = Counter()
    identified_counts = 0
    for hashval, lineages in hash_assign.items():
        if lineages:
            identified_counts += 1
            for lineage in lineages:
                lineage = utils.pop_to_rank(lineage, 'genus')
                counts[lineage] += 1

    # make sure it's strain or species level
    genome_lineage, count = next(iter(counts.most_common()))
    f_major = count / identified_counts
    print(f'{f_major*100:.1f}% of hashes identify as {pretty_print_lineage(genome_lineage)}')

    if args.lineage:
        provided_lin = args.lineage.split(';')
        provided_lin = [ LineagePair(rank, name) for (rank, name) in zip(sourmash.lca.taxlist(), provided_lin) ]
        print(f'provided lineage: {sourmash.lca.display_lineage(provided_lin)}')

        if utils.is_lineage_match(provided_lin, genome_lineage, 'genus'):
            print(f'XXX agree')
        else:
            print(f'XXX disagree')
            print('XXX', sourmash.lca.display_lineage(provided_lin))
            print('XXX', sourmash.lca.display_lineage(genome_lineage))

    if genome_lineage[-1].rank != 'genus':
        print(f'rank of major assignment is f{genome_lineage[-1].rank}; quitting')
        comment = f'rank of major assignment is f{genome_lineage[-1].rank}; needs to be genus'
        if args.summary:
            with open(args.summary, 'wt') as fp:
                w = csv.writer(fp)
                w.writerow([args.genome] + [""]*14 + [comment])
        open(args.report, 'wt').close()
        open(args.clean, 'wt').close()
        open(args.dirty, 'wt').close()
        sys.exit(0)

    # report everything...
    report_fp = open(args.report, 'wt')
    print(f'{f_major*100:.1f}% of hashes identify as {pretty_print_lineage(genome_lineage)}', file=report_fp)
    print(f'({identified_counts} identified hashes, {count} in most common)', file=report_fp)
    if f_major < 0.8:
        print(f'** WARNING ** majority lineage is less than 80% of assigned lineages. Beware!', file=report_fp)

    print(f'\n** hashval lineage counts for genome - {len(hash_assign)} => {len(hash_assign)*scaled/1000:.0f} kb', file=report_fp)
    for lin, count in counts.most_common():
        print(f'   {count*scaled/1000:.0f} kb {pretty_print_lineage(lin)}', file=report_fp)
        print(sourmash.lca.display_lineage(lin))
    print('', file=report_fp)

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

    print(f'pass 2: reading contigs from {args.genome}')
    print(f'**\n** walking through contigs:\n**\n', file=report_fp)
    for n, record in enumerate(screed.open(args.genome)):
        # make a new minhash and start examining it.
        mh = empty_mh.copy_and_clear()
        mh.add_sequence(record.sequence, force=True)

        clean = True               # default to clean
        if not mh:                 # no hashes?
            missed_n += 1
            missed_bp += len(record.sequence)

        if mh and len(mh) >= 2:
            clean = check_gather(record, mh, genome_lineage, lca_db, ldb, report_fp)
            if not clean:
                n_reason_3 += 1

        # did we find a dirty contig in step 1? if NOT, go into LCA style
        # approaches.
        if mh and clean:
            
            # get _all_ of the hash taxonomy assignments for this contig
            ctg_assign = gather_assignments(mh.get_mins(), None, [lca_db], ldb)

            ctg_tax_assign = count_lca_for_assignments(ctg_assign)
            if ctg_tax_assign:
                # get top assignment for contig.
                ctg_lin, lin_count = next(iter(ctg_tax_assign.most_common()))

                # assignment outside of genus? dirty!
                if ctg_lin[-1].rank not in ('species', 'strain', 'genus'):
                    clean = False
                    n_reason_1 += 1
                    print(f'\n---- contig {record.name} ({len(record.sequence)/1000:.0f} kb)', file=report_fp)
                    print(f'contig dirty, REASON 1 - contig LCA is above genus\nlca rank is {ctg_lin[-1].rank}',
                          file=report_fp)
                    print('', file=report_fp)
                elif not utils.is_lineage_match(genome_lineage, ctg_lin, 'genus'):
                    clean = False
                    n_reason_2 += 1
                    print('', file=report_fp)
                    print(f'---- contig {record.name} ({len(record.sequence)/1000:.0f} kb)', file=report_fp)
                    print(f'contig dirty, REASON 2 - contig lineage is not a match to genome\'s genus\nlineage is {pretty_print_lineage(ctg_lin)}',
                          file=report_fp)

                # summary reporting --
                if not clean:
                    report_lca_summary(report_fp, ctg_tax_assign,
                                       ctg_assign, scaled)

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
    print(f'kept {clean_n} contigs containing {int(clean_bp/1000)} kb.',
          file=report_fp)
    print(f'removed {dirty_n} contigs containing {int(dirty_bp/1000)} kb.',
          file=report_fp)
    print(f'{missed_n} contigs ({int(missed_bp/1000)} kb total) had no hashes, so counted as clean', file=report_fp)

    # look at what our database says about remaining contamination,
    # across all "clean" contigs. (Need to dig into this more to figure
    # out exactly why we still have any :)
    print(f'\nbreakdown of clean contigs w/gather:', file=report_fp)
    # CTB: add breakdown of dirty contigs?

    # report gather breakdown of clean signature
    clean_mh = clean_out.minhash
    first_match = do_gather_breakdown(clean_mh, lca_db, report_fp)

    # get genome size and match lineage of primary match
    nearest_size = 0
    match_lineage = ""
    if first_match:
        nearest_size = len(first_match.minhash) * first_match.minhash.scaled
        ident = get_ident(first_match)
        match_lineage = ldb.ident_to_lineage[ident]
        ratio = round(clean_bp / nearest_size, 2)

    # write out a one line summary?
    if args.summary:
        comment = ""

        with open(args.summary, 'wt') as fp:
            full_lineage = sourmash.lca.display_lineage(match_lineage)
            short_lineage = pretty_print_lineage(match_lineage)
            w = csv.writer(fp)
            w.writerow([args.genome, short_lineage, full_lineage,
                        nearest_size, ratio, clean_bp,
                        clean_n, dirty_n, dirty_bp,
                        missed_n, missed_bp, f_major,
                        n_reason_1, n_reason_2, n_reason_3,
                        comment])


if __name__ == '__main__':
    main()
