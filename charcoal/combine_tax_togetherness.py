#! /usr/bin/env python
# next steps -
#   load hashes into labels here, and add tax info
#   represent as tree, decoraet internal nodes
#   fun, profit.
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
from pickle import load, dump


def load_and_normalize(filename, delete_empty=False):
    mat = genfromtxt(filename, delimiter=',')
    assert mat.shape[0] == 347            # number of metagenomes
    n_hashes = mat.shape[1]
    n_orig_hashes = n_hashes

    # go through and normalize all the sample-presence vectors for each hash;
    # track those with all 0s for later removal.
    to_delete = []
    for i in range(n_hashes):
        if sum(mat[:, i]):
            mat[:, i] /= math.sqrt(np.dot(mat[:, i], mat[:, i]))
        else:
            to_delete.append(i)

    if delete_empty:
        # remove all columns with zeros
        print('removing {} null presence vectors'.format(len(to_delete)))
        for row_n in reversed(to_delete):
            mat = np.delete(mat, row_n, 1)

        assert mat.shape[1] == n_hashes - len(to_delete)

    n_hashes = mat.shape[1]

    # construct distance matrix using angular distance
    D = np.zeros((n_hashes, n_hashes))
    for i in range(n_hashes):
        for j in range(n_hashes):
            cos_sim = np.dot(mat[:, i], mat[:, j])
            cos_sim = min(cos_sim, 1.0)
            ang_sim = 1 - 2*math.acos(cos_sim) / math.pi
            D[i][j] = ang_sim

    # done!
    return D, n_orig_hashes


def do_cluster(mat, hashes_to_tax):
    n_hashes = mat.shape[1]
    assert len(hashes_to_tax) == n_hashes
    
    Y = sch.linkage(mat, method='complete')
    rootnode, nodelist = sch.to_tree(Y, rd=True)

    hashlist = sorted(hashes_to_tax)
    node_id_to_tax = {}
    for pos in range(n_hashes):
        node = nodelist[pos]
        assert node.get_id() == pos
        tax_info, reason = hashes_to_tax[hashlist[pos]]

        x = set()
        if tax_info:
            x.add(tax_info)
        node_id_to_tax[node.get_id()] = x

    def traverse_and_label(node, indent=' '):
        if node.is_leaf():
            tax_id = node_id_to_tax[node.get_id()]
            x = set()
            if tax_id:
                x.update(tax_id)
            return x
        else:
            l = traverse_and_label(node.get_left(), indent=indent + ' ')
            r = traverse_and_label(node.get_right(), indent=indent + ' ')

            x = l.union(r)
            node_id_to_tax[node.get_id()] = x

            return x

    traverse_and_label(rootnode)

    return rootnode, nodelist, node_id_to_tax

    def print_lca(node):
        node_id = node.get_id()
        tax_set = node_id_to_tax[node_id]

        if not tax_set:
            return "-nada-", "none"
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
        lca_str = lca_str.replace(';', '+').replace(' ', '_').replace('_', '')
        is_leaf = ' '
        if node.is_leaf():
            is_leaf = '*'
        print('YYY', indent, node.get_id(), is_leaf, lca_str, reason)
        if node.is_leaf():
            return '{}'.format(lca_str)
        else:
            l = traverse(node.get_left(), indent=indent + ' ')
            r = traverse(node.get_right(), indent=indent + ' ')
            return "({},{}){}".format(l, r, lca_str)

    newick_tree = traverse(rootnode) + ';'

    # we should have 2*n_hashes - 1 clusters
    assert len(nodelist) == 2*n_hashes - 1

    # ok. so this linkage matrix has entries where each cluster k
    # is composed of two IDs, Y[k][0] and Y[k][1]. If an ID has
    # value (0..n) where n is original observation, then it is an
    # original observation. Otherwise it is a non-leaf node ID.
    
    # or, in reverse:
    # every entry i in Y[*,0] and Y[*,1] refers to an entry in nodelist.

    nodelist_to_linkage = [-1]*len(nodelist)
    for k in range(n_hashes - 1):
        id_1 = int(Y[k][0])
        id_2 = int(Y[k][1])

        assert nodelist_to_linkage[id_1] == -1
        assert nodelist_to_linkage[id_2] == -1
        nodelist_to_linkage[id_1] = k
        nodelist_to_linkage[id_2] = k

    # the last node in the list is the only one without any original
    # observations belonging to it.
    last_node = nodelist[-1]
    assert last_node.get_left().get_id() >= n_hashes
    assert last_node.get_right().get_id() >= n_hashes

#    Z = sch.dendrogram(Y, orientation='left')
#    print(Z['leaves'])
#    print(Z['ivl'])

    return newick_tree


def main():
    p = argparse.ArgumentParser()
    p.add_argument('matrix_csv')
    p.add_argument('taxhashes')
    p.add_argument('--newick', default=None)
    p.add_argument('--pickle-tree', default=None)
    args = p.parse_args()

    mat, n_orig_hashes = load_and_normalize(args.matrix_csv)

    with open(args.taxhashes, 'rb') as fp:
        hashes_to_tax = load(fp)

    for k in hashes_to_tax:
        lca, reason = hashes_to_tax[k]

    assert n_orig_hashes == len(hashes_to_tax)
    assert mat.shape[1] == len(hashes_to_tax)

    rootnode, nodelist, node_id_to_tax = do_cluster(mat, hashes_to_tax)

    if args.pickle_tree:
        with open(args.pickle_tree, 'wb') as fp:
            dump((rootnode, nodelist, node_id_to_tax), fp)

#    newick_tree = do_cluster(mat, hashes_to_tax)

#    if args.newick:
#        with open(args.newick, 'wt') as fp:
#            fp.write(newick_tree)


if __name__ == '__main__':
    main()
