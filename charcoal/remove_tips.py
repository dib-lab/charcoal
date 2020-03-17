#! /usr/bin/env python
"""
Remove bad tips based on taxonomy.

A zeroth attempt to just get something basic working for contaminant removal.

This is mostly just internal error checking...
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
    p.add_argument('pickle_tree', help='output of combine_tax_togetherness')
    p.add_argument('hashes', help='output of process_genome')
    p.add_argument('--rm-hashes', help='output hashes to remove')
    args = p.parse_args()

    with open(args.pickle_tree, 'rb') as fp:
        (rootnode, nodelist, node_id_to_tax) = load(fp)

    # find majority across leaves
    leaf_tax = collections.Counter()
    for node in nodelist:
        if node.is_leaf():
            tax_set = node_id_to_tax[node.get_id()]
            assert len(tax_set) <= 1
            if tax_set:
                p = list(tax_set)[0]
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

    rm_leaves = set()
    for node in nodelist:
        if node.is_leaf():
            tax_set = node_id_to_tax[node.get_id()]
            assert len(tax_set) <= 1
            if tax_set:
                node_lineage = list(tax_set)[0]
                if not is_lineage_match(node_lineage, most_common, 'order'):
                    print(lca_utils.display_lineage(node_lineage))
                    rm_leaves.add(node.get_id())
    print(rm_leaves)

    with open(args.hashes, 'rb') as fp:
        hash_to_lengths = load(fp)

    rm_hashes = set()
    for hashpos, hash in enumerate(hash_to_lengths):
        if hashpos in rm_leaves:
            rm_hashes.add(hash)
            rm_leaves.remove(hashpos)
    assert not rm_leaves

    with open(args.rm_hashes, 'wt') as fp:
        print("\n".join([str(h) for h in sorted(rm_hashes) ]), file=fp)
        

if __name__ == '__main__':
    main()
