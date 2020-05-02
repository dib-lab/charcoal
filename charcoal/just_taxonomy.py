#! /usr/bin/env python
"""
Remove bad contigs based solely on taxonomy.

CTB TODO:
* optionally eliminate contigs with no taxonomy
"""
import argparse
import gzip
from collections import Counter

import screed

import sourmash
from sourmash.lca.command_index import load_taxonomy_assignments
from sourmash.lca import LCA_Database

from . import utils


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


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--genome', help='genome file', required=True)
    p.add_argument('--lineages_csv', help='lineage spreadsheet', required=True)
    p.add_argument('--matches_sig', help='all relevant matches', required=True)
    p.add_argument('--clean', help='cleaned contigs', required=True)
    p.add_argument('--dirty', help='dirty contigs', required=True)
    p.add_argument('--report', help='report output', required=True)
    args = p.parse_args()

    tax_assign, _ = load_taxonomy_assignments(args.lineages_csv,
                                              start_column=3)
    print(f'loaded {len(tax_assign)} tax assignments.')

    with open(args.matches_sig, 'rt') as fp:
        siglist = list(sourmash.load_signatures(fp))

    if not siglist:
        print('no matches for this genome, exiting.')
        sys.exit(-1)

    empty_mh = siglist[0].minhash.copy_and_clear()
    ksize = empty_mh.ksize
    scaled = empty_mh.scaled

    lca_db = LCA_Database(ksize=ksize, scaled=scaled)

    for ss in siglist:
        ident = ss.name()
        ident = ident.split()[0]
        ident = ident.split('.')[0]
        lineage = tax_assign[ident]

        lca_db.insert(ss, ident=ident, lineage=lineage)

    print(f'loaded {len(siglist)} signatures & created LCA Database')

    print(f'pass 1: reading contigs from {args.genome}')
    entire_mh = empty_mh.copy_and_clear()
    for n, record in enumerate(screed.open(args.genome)):
        entire_mh.add_sequence(record.sequence, force=True)

    # get all of the hash taxonomy assignments for this contig
    hash_assign = sourmash.lca.gather_assignments(entire_mh.get_mins(),
                                                  [lca_db])

    # count them and find major
    counts = Counter()
    identified_counts = 0
    for hashval, lineages in hash_assign.items():
        if lineages:
            identified_counts += 1
            for lineage in lineages:
                counts[lineage] += 1

    # make sure it's strain or species level
    assign, count = next(iter(counts.most_common()))
    f_major = count / len(hash_assign)
    print(f'{f_major*100:.1f}% of hashes identify as {pretty_print_lineage(assign)}')
    if assign[-1].rank not in ('species', 'strain'):
        print(f'rank of major assignment is f{assign[-1].rank}; quitting')
        sys.exit(-1)

    report_fp = open(args.report, 'wt')
    print(f'{f_major*100:.1f}% of hashes identify as {pretty_print_lineage(assign)}', file=report_fp)
    print(f'({identified_counts} identified hashes, {count} in most common)', file=report_fp)

    clean_fp = gzip.open(args.clean, 'wt')
    dirty_fp = gzip.open(args.dirty, 'wt')

    # now, find disagreeing contigs.
    dirty_bp = clean_bp = 0
    dirty_n = clean_n = 0
    
    print(f'pass 2: reading contigs from {args.genome}')
    for n, record in enumerate(screed.open(args.genome)):
        clean = True               # default to clean

        mh = empty_mh.copy_and_clear()
        mh.add_sequence(record.sequence, force=True)

        if mh:
            # get all of the hash taxonomy assignments for this contig
            ctg_assign = sourmash.lca.gather_assignments(mh.get_mins(),
                                                         [lca_db])

            ctg_tax_assign = sourmash.lca.count_lca_for_assignments(ctg_assign)
            if ctg_tax_assign:
                ctg_lin, lin_count = next(iter(ctg_tax_assign.most_common()))

                # assignment outside of genus? dirty!
                if ctg_lin[-1].rank not in ('species', 'strain', 'genus'):
                    clean = False
                    print(f'---- contig {record.name}')
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
        else:
            dirty_fp.write(f'>{record.name}\n{record.sequence}\n')
            dirty_n += 1
            dirty_bp += len(record.sequence)

    print('--------------', file=report_fp)
    print(f'kept {clean_n} contigs containing {int(clean_bp)/1000} kb.',
          file=report_fp)
    print(f'removed {dirty_n} contigs containing {int(dirty_bp)/1000} kb.',
          file=report_fp)


if __name__ == '__main__':
    main()
