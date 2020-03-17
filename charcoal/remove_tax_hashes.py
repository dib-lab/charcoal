#! /usr/bin/env python
"""
Remove bad hashes based solely on taxonomy.

A test to get something basic working for contaminant removal.
"""
import argparse
from sourmash.lca import lca_utils

import collections
from pickle import load
import pprint


def is_lineage_match(lin_a, lin_b, rank):
    for a, b in zip(lin_a, lin_b):
        assert a.rank == b.rank
        if a.rank == rank:
            if a == b:
                return 1
        if a != b:
            return 0

    return 0


def main():
    p = argparse.ArgumentParser()
    p.add_argument('tax_hashes', help='output of genome_shred_to_tax')
    p.add_argument('--rm-hashes', help='output hashes to remove')
    args = p.parse_args()

    with open(args.tax_hashes, 'rb') as fp:
        hashes_to_tax = load(fp)

    # find majority across leaves
    leaf_tax = collections.Counter()
    for hashval, (tax, reason) in hashes_to_tax.items():
        if tax:
            p = tax
            p = list(p)
            while p and p[-1].rank != 'order':
                p.pop()

            if p and p[-1].rank == 'order':
                leaf_tax[tuple(p)] += 1

    for k, v in leaf_tax.most_common():
        print('lineage {} has count {}'.format(lca_utils.display_lineage(k), v))
    print('')

    most_common, most_common_count = next(iter(leaf_tax.most_common(1)))
    print('removing all but {}'.format(lca_utils.display_lineage(most_common)))

    rm_hashes = set()
    for hashval, (tax, reason) in hashes_to_tax.items():
        if tax:
            if not is_lineage_match(tax, most_common, 'order'):
                print(lca_utils.display_lineage(tax))
                rm_hashes.add(hashval)
    print(rm_hashes)

    with open(args.rm_hashes, 'wt') as fp:
        print("\n".join([str(h) for h in rm_hashes ]), file=fp)
        

if __name__ == '__main__':
    main()
