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

    def find_lca(node_id):
        tax_set = node_id_to_tax[node_id]
        if not tax_set:
            return None

        tree = lca_utils.build_tree(tax_set)
        lca, n_children = lca_utils.find_lca(tree)
        return lca

    rm_nodes = set()
    def get_lca_str(node_id):
        tax_set = node_id_to_tax[node_id]
        if not tax_set:
            return "- none -"
        tree = lca_utils.build_tree(tax_set)
        lca, n_children = lca_utils.find_lca(tree)

        if lca and lca[-1].rank in ('superkingdom', 'phylum', 'class'):
            if len(tax_set) > 1:
                print('-'*20)
                lca_str = lca_utils.display_lineage(lca, include_strain=False)
                print('LCA for node {}; {} children. {}'.format(node_id, len(tax_set), lca_str))
                for k in tax_set:
                    print('   * ', lca_utils.display_lineage(k, include_strain=False, truncate_empty=False))

            node = nodelist[node_id]
            if not node.is_leaf():
                l = node.get_left()
                l_lca = find_lca(l.get_id())
                r = node.get_right()
                r_lca = find_lca(r.get_id())
                keep_going = 1
                if l_lca and l_lca[-1].rank in ('superkingdom', 'phylum', 'class'):
                    keep_going = 0
                if r_lca and r_lca[-1].rank in ('superkingdom', 'phylum', 'class'):
                    keep_going = 0

                if keep_going:
                    print('xxx 1', lca_utils.display_lineage(l_lca, include_strain=False, truncate_empty=False))
                    print('xxx 2', lca_utils.display_lineage(r_lca, include_strain=False, truncate_empty=False))

                    rm = 0
                    rm_node_id = None
                    if is_lineage_match(l_lca, most_common, 'order'):
                        print('(xxx keeping first)')
                        rm = 1
                        rm_node_id = l.get_id()
                    if is_lineage_match(r_lca, most_common, 'order'):
                        print('(xxx keeping second)')
                        assert not rm
                        rm = 1
                        rm_node_id = r.get_id()
                    assert rm
                    rm_nodes.add(rm_node_id)
                

        lca = list(lca)

        # find species, or next best thing
        for i in range(len(lca)):
            if lca[-1] and lca[-1].rank == 'strain':
                lca.pop()
            elif lca[-1] and lca[-1].rank == 'species':
                lca_str = lca[-1].name
                break
            elif lca[-1]:
                lca_str = "{}={}".format(lca[-1].rank, lca[-1].name)
                break

        return lca_str

    lca_str_set = set()
    for k in node_id_to_tax:
        lca_str = get_lca_str(k)
        lca_str_set.add(lca_str)

    taxon_namespace = dendropy.TaxonNamespace(lca_str_set)
    tree = dendropy.Tree(taxon_namespace=taxon_namespace)
    
    def traverse(node):
        lca_str = get_lca_str(node.get_id())
        ch = dendropy.Node(edge_length=1)
        ch.taxon = taxon_namespace.get_taxon(lca_str)
        
        if node.is_leaf():
            return ch
        else:
            l = traverse(node.get_left())
            r = traverse(node.get_right())

            ch.set_child_nodes([l, r])
            return ch

    tree.seed_node.add_child(traverse(rootnode))

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
    for hashpos, hash in enumerate(hash_to_lengths):
        if hashpos in rm_leaves:
            print(hashpos, hash)
            rm_hashes.add(hash)
            rm_leaves.remove(hashpos)
    assert not rm_leaves

    with open(args.rm_hashes, 'wt') as fp:
        print("\n".join([str(h) for h in rm_hashes ]), file=fp)
        
#    with open(args.newick_out, 'wt') as fp:
#        print(tree.as_string("newick"), file=fp)

#    print(tree.as_ascii_plot())


if __name__ == '__main__':
    main()
