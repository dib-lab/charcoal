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

from utils import is_lineage_match, pop_to_rank


def main():
    p = argparse.ArgumentParser()
    p.add_argument('tax_hashes', help='output of genome_shred_to_tax')
    p.add_argument('--rm-hashes', help='output hashes to remove')
    args = p.parse_args()

    with open(args.tax_hashes, 'rb') as fp:
        hashes_to_tax = load(fp)

    # find majority across leaves
    leaf_tax = collections.Counter()
    for hashval, tax in hashes_to_tax.items():
        if tax:
            p = pop_to_rank(tax, 'order')
            if p:
                leaf_tax[tuple(p)] += 1

    for k, v in leaf_tax.most_common():
        print('lineage {} has count {}'.format(lca_utils.display_lineage(k), v))
    print('')

    most_common, most_common_count = next(iter(leaf_tax.most_common(1)))
    print('removing all but {}'.format(lca_utils.display_lineage(most_common)))

    # find all hashes belonging to "bad" (all but most common) tax.
    rm_hashes = set()
    for hashval, tax in hashes_to_tax.items():
        if tax:
            if not is_lineage_match(tax, most_common, 'order'):
                print(lca_utils.display_lineage(tax))
                rm_hashes.add(hashval)
    print(rm_hashes)

    # save!
    with open(args.rm_hashes, 'wt') as fp:
        print("\n".join([str(h) for h in sorted(rm_hashes) ]), file=fp)
        

if __name__ == '__main__':
    main()
