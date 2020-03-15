#! /usr/bin/env python
"""
Write out a Newick tree, using dendropy, from output of 'combine_tax_togetherness'
"""
import argparse
from sourmash.lca import lca_utils

import collections
from pickle import load
import dendropy

def main():
    p = argparse.ArgumentParser()
    p.add_argument('pickle_tree', help='output of combine_tax_togetherness')
    p.add_argument('newick_out', help='output newick tree')
    args = p.parse_args()

    with open(args.pickle_tree, 'rb') as fp:
        (rootnode, nodelist, node_id_to_tax) = load(fp)

    def get_lca_str(node_id):
        tax_set = node_id_to_tax[node_id]
        if not tax_set:
            return "- none -"
        tree = lca_utils.build_tree(tax_set)
        lca, reason = lca_utils.find_lca(tree)

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
    
    def traverse(node, indent=' '):
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
    with open(args.newick_out, 'wt') as fp:
        print(tree.as_string("newick"), file=fp)

    print(tree.as_ascii_plot())


if __name__ == '__main__':
    main()
