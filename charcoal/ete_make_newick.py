#! /usr/bin/env python
"""
Write out a Newick tree, using ete3, for output of 'combine-tax-togetherness.py'
"""
import argparse
import numpy as np
from numpy import genfromtxt
import math
from scipy.cluster.hierarchy import dendrogram, linkage
import pprint
from sourmash.lca import lca_utils

from matplotlib import pyplot as plt
import pylab
import scipy.cluster.hierarchy as sch
import collections
from pickle import load
from ete3 import Tree

def main():
    p = argparse.ArgumentParser()
    p.add_argument('pickle_tree')
    p.add_argument('newick_out')
    args = p.parse_args()

    with open(args.pickle_tree, 'rb') as fp:
        (rootnode, nodelist, node_id_to_tax) = load(fp)

    def print_lca(node):
        node_id = node.get_id()
        tax_set = node_id_to_tax[node_id]

        if not tax_set:
            return "unknown", "none"
        tree = lca_utils.build_tree(tax_set)
        lca, reason = lca_utils.find_lca(tree)

        lca = list(lca)

        # find species, or next best thing
        for i in range(len(lca)):
            if lca[-1] and lca[-1].rank == 'strain':
                lca.pop()
            elif lca[-1] and lca[-1].rank == 'species':
                return lca[-1].name, reason
            elif lca[-1]:
                return "{}={}".format(lca[-1].rank, lca[-1].name), reason

        return lca_utils.display_lineage(lca, truncate_empty=True), reason

    def traverse(node, indent=' '):
        lca_str, reason = print_lca(node)
        lca_str = lca_str.replace(';', '+')# .replace(' ', '_').replace('_', '')
        is_leaf = ' '
        if node.is_leaf():
            is_leaf = '*'
        print('YYY', indent, node.get_id(), is_leaf, lca_str, reason)
        if node.is_leaf():
            tn = Tree()
            tn.add_child(name=lca_str)
            return tn
        else:
            tn = Tree()
            B = tn.add_child(name=lca_str)
            B.add_child(traverse(node.get_left(), indent=indent + ' '))
            B.add_child(traverse(node.get_right(), indent=indent + ' '))
            return B

    T = traverse(rootnode)
    with open(args.newick_out, 'wt') as fp:
        print(T.write(format=0), file=fp)


if __name__ == '__main__':
    main()
