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
    p.add_argument('--tax-hashes')
    p.add_argument('-o', '--output')
    args = p.parse_args()

    assert args.tax_hashes
    assert args.output

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


if __name__ == '__main__':
    sys.exit(main())
