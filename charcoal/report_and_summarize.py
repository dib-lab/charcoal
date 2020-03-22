#! /usr/bin/env python
"""
Report a summary of what charcoal has seen and is doing.
"""

import sys
import argparse
import os
from pickle import load
from collections import Counter

import screed
import sourmash
from sourmash.lca import lca_utils

import utils                              # charcoal utils


def main():
    p = argparse.ArgumentParser()
    p.add_argument('genome')
    p.add_argument('--tax-hashes', help='output of genome_shred_to_tax')
    p.add_argument('--matrix', help='output of match_metagenomes')
    p.add_argument('-o', '--output')
    args = p.parse_args()

    assert args.tax_hashes
    assert args.output
    assert args.matrix

    n_contigs = 0
    sum_bp = 0
    for record in screed.open(args.genome):
        n_contigs += 1
        sum_bp += len(record.sequence)

    sum_mbp = sum_bp / 1e6

    outfp = open(args.output, 'wt')
    print(f"""\
# Report: {os.path.basename(args.genome)}

Genome file: {args.genome}

{sum_bp / 1e6:.1f} Mbp in {n_contigs} contigs.
""", file=outfp)

    ####

    with open(args.tax_hashes, 'rb') as fp:
        hashes_to_tax = load(fp)

    total_contigs = 0
    sum_bp = 0
    missed_contigs = 0
    sum_missed_bp = 0
    mh = sourmash.MinHash(n=0, ksize=hashes_to_tax.ksize,
                          scaled=hashes_to_tax.scaled)
    shredder = utils.GenomeShredder(args.genome, hashes_to_tax.fragment_size)
    for name, seq, start, end in shredder:
        total_contigs += 1
        sum_bp += len(seq)
        
        mh.add_sequence(seq, force=True)
        if not mh:
            sum_missed_bp += len(seq)
            missed_contigs += 1

    #assert total_contigs == len(hashes_to_tax) @CTB WTF

    print(f"""\
## Contigs & fragments report

Fragment size: {hashes_to_tax.fragment_size}
Number of fragments recorded w/hashes: {len(hashes_to_tax)} 
sourmash ksize: {hashes_to_tax.ksize}
sourmash scaled: {hashes_to_tax.scaled}

Hashing missed {sum_missed_bp/1000:.2f} kbp in {missed_contigs} contigs.
""", file=outfp)

    ### order

    lca_count_order = Counter()
    for v in hashes_to_tax.d.values():
        v2 = utils.pop_to_rank(v, 'order')
        v = tuple(v2)
        lca_count_order[v2] += 1

    print(f"""\
## Taxonomy report:

LCA database: {hashes_to_tax.lca_db_file}

Order or above: {len(lca_count_order)} distinct taxonomic matches in the genome fragments.
""", file=outfp)

    for lca, cnt in lca_count_order.most_common():
        pcnt = cnt / len(hashes_to_tax) * 100.
        print(f"""\
* {cnt} fragments ({pcnt:.1f}%), {lca_utils.display_lineage(lca, truncate_empty=False)}""",
              file=outfp)


    lca_count_genus = Counter()
    for v in hashes_to_tax.d.values():
        v2 = utils.pop_to_rank(v, 'genus')
        v = tuple(v2)
        lca_count_genus[v2] += 1

    print(f"""\

----

Genus or above: {len(lca_count_genus)} distinct taxonomic matches in the genome fragments.
""", file=outfp)

    for lca, cnt in lca_count_genus.most_common():
        pcnt = cnt / len(hashes_to_tax) * 100.
        print(f"""\
* {cnt} fragments ({pcnt:.1f}%), {lca_utils.display_lineage(lca, truncate_empty=False)}""",
              file=outfp)

    rank_count = Counter()
    for lca, cnt in lca_count_genus.most_common():
        rank = '(none)'
        if lca:
            rank = lca[-1].rank
        rank_count[rank] += 1

    print(f"\nGenus or above: {len(rank_count)} ranks in tax classifications.", file=outfp)
    for rank, cnt in rank_count.most_common():
        print(f"* {cnt} at rank '{rank}'", file=outfp)

    ###

    print('calculating distance matrix from', args.matrix)
    with open(args.matrix, 'rb') as fp:
        matrix_obj = load(fp)

    matrix = matrix_obj.mat

    n_metagenomes, n_hashes = matrix.shape

    n_not_present = 0
    for i in range(n_metagenomes):
        if not sum(matrix[i, :]):
            n_not_present += 1

    n_empty_hashes = 0
    for j in range(n_hashes):
        if not sum(matrix[:, j]):
            n_empty_hashes += 1

    # @CTB include plot?
    print(f"""
## Metagenome togetherness

Across {n_metagenomes} metagenomes, {n_hashes - n_empty_hashes} of
{n_hashes} hashes are present in at least one sample.

Of {n_metagenomes} metagenomes, this genome has no overlap with
{n_not_present} of them.

""", file=outfp)


if __name__ == '__main__':
    sys.exit(main())
