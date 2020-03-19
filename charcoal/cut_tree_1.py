#! /usr/bin/env python
"""
Cut togetherness+tax tree at transitions genus/species/strain to order.

A first attempt to just get something basic working for contaminant removal. :)
"""
import argparse
from sourmash.lca import lca_utils

import collections
from pickle import load
import pprint

import dendropy

from utils import is_lineage_match, pop_to_rank


def make_lca(node, node_id_to_tax):
    tax_set = node_id_to_tax[node.get_id()]
    if not tax_set:
        return None

    tree = lca_utils.build_tree(tax_set)
    lca, n_children = lca_utils.find_lca(tree)
    return lca


def query_cut_node(node, node_id_to_tax, most_common):
    # should we cut this node?
    lca = make_lca(node, node_id_to_tax)

    # no lca? ignore.
    if not lca:
        return False

    # high level node? ignore. confusion below.
    if lca[-1].rank in ('superkingdom', 'phylum', 'class'):
        return False

    # matches most common lineage? do not cut.
    if is_lineage_match(lca, most_common, 'order'):
        return False

    # cut!
    return True


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
                p = pop_to_rank(list(tax_set)[0], 'order')

                if p and p[-1].rank == 'order':
                    leaf_tax[tuple(p)] += 1

    for k, v in leaf_tax.most_common():
        print('lineage {} has count {}'.format(lca_utils.display_lineage(k), v))
    print('')

    most_common, most_common_count = next(iter(leaf_tax.most_common(1)))
    print('removing all but {}'.format(lca_utils.display_lineage(most_common)))

    rm_nodes = set()
    for node in nodelist:
        if query_cut_node(node, node_id_to_tax, most_common):
            rm_nodes.add(node.get_id())
    print(rm_nodes)

    rm_leaves = set()
    def traverse_add_rm(node):
        if node.is_leaf():
            rm_leaves.add(node.get_id())
        else:
            l = traverse_add_rm(node.get_left())
            r = traverse_add_rm(node.get_right())

    for node_id in rm_nodes:
        traverse_add_rm(nodelist[node_id])

    print(rm_leaves)

    with open(args.hashes, 'rb') as fp:
        hash_to_lengths = load(fp)

    rm_hashes = set()
    for hashpos, hash in enumerate(sorted(hash_to_lengths)):
        if hashpos in rm_leaves:
            rm_hashes.add(hash)
            rm_leaves.remove(hashpos)
    assert not rm_leaves

    with open(args.rm_hashes, 'wt') as fp:
        print("\n".join([str(h) for h in sorted(rm_hashes) ]), file=fp)
        

if __name__ == '__main__':
    main()
