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
from sourmash.lca import LCA_Database

from . import utils
from . import lineage_db
from .lineage_db import LineageDB


def get_idents_for_hashval(lca_db, hashval):
    idx_list = lca_db.hashval_to_idx.get(hashval, [])
    for idx in idx_list:
        ident = lca_db.idx_to_ident[idx]
        yield ident


def gather_assignments(hashvals, rank, dblist, ldb):
    """
    Gather assignments from across all the databases for all the hashvals.
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
    ident = sig.name()
    ident = ident.split()[0]
    ident = ident.split('.')[0]
    return ident
    

def main():
    p = argparse.ArgumentParser()
    p.add_argument('--genome', help='genome file', required=True)
    p.add_argument('--lineages_csv', help='lineage spreadsheet', required=True)
    p.add_argument('--matches_sig', help='all relevant matches', required=True)
    p.add_argument('--clean', help='cleaned contigs', required=True)
    p.add_argument('--dirty', help='dirty contigs', required=True)
    p.add_argument('--report', help='report output', required=True)
    p.add_argument('--summary', help='CSV one line output')
    args = p.parse_args()

    tax_assign, _ = load_taxonomy_assignments(args.lineages_csv,
                                              start_column=3)
    print(f'loaded {len(tax_assign)} tax assignments.')

    with open(args.matches_sig, 'rt') as fp:
        siglist = list(sourmash.load_signatures(fp))

    if not siglist:
        print('no matches for this genome, exiting.')
        sys.exit(-1)

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
        lineage = tax_assign[ident]

        lca_db.insert(ss, ident=ident)
        ldb.insert(ident, lineage)

    print(f'loaded {len(siglist)} signatures & created LCA Database')

    print(f'pass 1: reading contigs from {args.genome}')
    entire_mh = empty_mh.copy_and_clear()
    for n, record in enumerate(screed.open(args.genome)):
        entire_mh.add_sequence(record.sequence, force=True)

    # get all of the hash taxonomy assignments for this contig
    hash_assign = gather_assignments(entire_mh.get_mins(), 'genus', [lca_db], ldb)

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
    assign, count = next(iter(counts.most_common()))
    f_major = count / identified_counts
    print(f'{f_major*100:.1f}% of hashes identify as {pretty_print_lineage(assign)}')
    if assign[-1].rank not in ('species', 'strain', 'genus'):
        print(f'rank of major assignment is f{assign[-1].rank}; quitting')
        sys.exit(-1)

    report_fp = open(args.report, 'wt')
    print(f'{f_major*100:.1f}% of hashes identify as {pretty_print_lineage(assign)}', file=report_fp)
    print(f'({identified_counts} identified hashes, {count} in most common)', file=report_fp)
    if f_major < 0.8:
        print(f'** WARNING ** majority lineage is less than 80% of assigned lineages. Beware!', file=report_fp)

    print(f'\n** hashval lineage counts for genome - {len(hash_assign)}', file=report_fp)
    for lin, count in counts.most_common():
        print(f'   {count} {pretty_print_lineage(lin)}', file=report_fp)

    clean_fp = gzip.open(args.clean, 'wt')
    dirty_fp = gzip.open(args.dirty, 'wt')

    # now, find disagreeing contigs.
    dirty_bp = clean_bp = 0
    dirty_n = clean_n = 0
    clean_mh = empty_mh.copy_and_clear()
    
    print(f'pass 2: reading contigs from {args.genome}')
    for n, record in enumerate(screed.open(args.genome)):
        clean = True               # default to clean

        mh = empty_mh.copy_and_clear()
        mh.add_sequence(record.sequence, force=True)

        if mh:
            # first, if there is more than one hash, use gather.
            if len(mh) >= 2:
                threshold_bp = mh.scaled*2
                results = lca_db.gather(sourmash.SourmashSignature(mh))
                if results:
                    match = results[0][1]
                    match_ident = get_ident(match)
                    lineage = tax_assign[match_ident]
                    if not utils.is_lineage_match(assign, lineage, 'genus'):
                        clean=False
                        print(f'contig {record.name} dirty, REASON 3\n   gather yields {pretty_print_lineage(lineage)}',
                              file=report_fp)
                        print('', file=report_fp)
                        print(f'---- contig {record.name}', file=report_fp)

        if mh and clean:
            
            # get all of the hash taxonomy assignments for this contig
            ctg_assign = gather_assignments(mh.get_mins(), None, [lca_db], ldb)

            ctg_tax_assign = count_lca_for_assignments(ctg_assign)
            if ctg_tax_assign:
                ctg_lin, lin_count = next(iter(ctg_tax_assign.most_common()))

                # assignment outside of genus? dirty!
                if ctg_lin[-1].rank not in ('species', 'strain', 'genus'):
                    clean = False
                    print(f'\n---- contig {record.name}')
                    print(f'dirty! {ctg_lin[-1].rank}')
                    print(f'contig {record.name} dirty, REASON 1\nrank is {ctg_lin[-1].rank}',
                          file=report_fp)
                    print('', file=report_fp)
                    print(f'---- contig {record.name}', file=report_fp)
                elif not utils.is_lineage_match(assign, ctg_lin, 'genus'):
                    clean = False
                    print(f'dirty! {ctg_lin}')
                    print('', file=report_fp)
                    print(f'---- contig {record.name}', file=report_fp)
                    print(f'contig {record.name} dirty, REASON 2\nlineage is {pretty_print_lineage(ctg_lin)}',
                          file=report_fp)

                # summary reporting

                if not clean:
                    ctg_counts = Counter()
                    for hashval, lineages in ctg_assign.items():
                        for lineage in lineages:
                            ctg_counts[lineage] += 1


                    print(f'\n** hashval lca counts', file=report_fp)
                    for lin, count in ctg_tax_assign.most_common():
                        print(f'   {count} {pretty_print_lineage(lin)}', file=report_fp)
                    print(f'\n** hashval lineage counts - {len(ctg_assign)}', file=report_fp)
                    for lin, count in ctg_counts.most_common():
                        print(f'   {count} {pretty_print_lineage(lin)}', file=report_fp)

        if clean:
            clean_fp.write(f'>{record.name}\n{record.sequence}\n')
            clean_n += 1
            clean_bp += len(record.sequence)

            clean_mh.add_sequence(record.sequence, force=True)
        else:
            dirty_fp.write(f'>{record.name}\n{record.sequence}\n')
            dirty_n += 1
            dirty_bp += len(record.sequence)

    print('--------------', file=report_fp)
    print(f'kept {clean_n} contigs containing {int(clean_bp)/1000} kb.',
          file=report_fp)
    print(f'removed {dirty_n} contigs containing {int(dirty_bp)/1000} kb.',
          file=report_fp)

    print(f'\nbreakdown of clean contigs w/gather:', file=report_fp)

    clean_sig = sourmash.SourmashSignature(clean_mh)
    while 1:
        results = lca_db.gather(clean_sig, threshold_bp=0)
        if not results:
            break

        (match, match_sig, _) = results[0]
        print(f'  {match*100:.3f}% - to {match_sig.name()}', file=report_fp)
        clean_mh.remove_many(match_sig.minhash.get_mins())
        clean_sig = sourmash.SourmashSignature(clean_mh)

    if args.summary:
        with open(args.summary, 'wt') as fp:
            w = csv.writer(fp)
            w.writerow([clean_n, clean_bp, dirty_n, dirty_bp,f_major])

if __name__ == '__main__':
    main()
