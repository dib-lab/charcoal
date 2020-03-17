#! /usr/bin/env python
"""
Combine taxonomy information and "togetherness" into one mishmash of data
structures.
"""
import argparse
import pprint
import scipy.cluster.hierarchy as sch
from pickle import load, dump

from sourmash.lca import lca_utils

from utils import load_and_normalize


def do_cluster(mat, hashes_to_tax):
    """
    Use scipy.cluster.hierarchy to cluster the matrix.
    """
    n_hashes = mat.shape[1]
    assert len(hashes_to_tax) == n_hashes

    # do the clustering...
    Y = sch.linkage(mat, method='complete')
    rootnode, nodelist = sch.to_tree(Y, rd=True)

    node_id_to_original = []              # len n_hashes
    for i in range(n_hashes - 1):
        if Y[i][0] < n_hashes:
            node_id_to_original.append(int(Y[i][0]))
        if Y[i][1] < n_hashes:
            node_id_to_original.append(int(Y[i][1]))

    print(len(node_id_to_original))
    print(node_id_to_original)

    # now, track the taxonomy <-> hashval <-> cluster node.
    hashlist = sorted(hashes_to_tax)
    node_id_to_tax = {}
    leaf_id_to_hash = {}

    # here we are assuming that the nodelist is in the same order as the
    # original labels. Is that right? NO, it is not.
    for pos in range(n_hashes):
        node = nodelist[pos]
        assert node.get_id() == pos
        original_id = node_id_to_original[pos]
        leaf_id_to_hash[pos] = hashlist[original_id]
        tax_info, reason = hashes_to_tax[hashlist[original_id]]

        x = set()
        if tax_info:
            x.add(tax_info)
        node_id_to_tax[pos] = x

    ###

    def traverse_and_label(node):
        "do recursive traversal from node and label with taxonomy LCA"
        if node.is_leaf():
            tax_id = node_id_to_tax[node.get_id()]
            x = set()
            if tax_id:
                x.update(tax_id)
            return x
        else:
            l = traverse_and_label(node.get_left())
            r = traverse_and_label(node.get_right())

            # taxonomy of a node is union of taxonomies below node.
            x = l.union(r)
            node_id_to_tax[node.get_id()] = x

            return x

    traverse_and_label(rootnode)

    ###

    return rootnode, nodelist, node_id_to_tax

def leftovers():                          # ignore me
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

    return newick_tree


def main():
    p = argparse.ArgumentParser()
    p.add_argument('matrix_csv', help='output of match_metagenomes')
    p.add_argument('hashes', help='output of process_genome')
    p.add_argument('taxhashes', help='output of genome_shred_to_tax')
    p.add_argument('--pickle-tree', default=None)
    args = p.parse_args()

    # output of match_metagenomes
    print('calculating distance matrix from', args.matrix_csv)
    mat, n_orig_hashes = load_and_normalize(args.matrix_csv)

    # output of process_genome
    with open(args.hashes, 'rb') as fp:
        hashes_to_length = load(fp)

    # output of genome_shred_to_tax
    with open(args.taxhashes, 'rb') as fp:
        hashes_to_tax = load(fp)

    # some basic validation
    assert n_orig_hashes == len(hashes_to_tax), "mismatch! was same --scaled used to compute these?"
    assert mat.shape[1] == len(hashes_to_tax)
    assert set(hashes_to_length) == set(hashes_to_tax)

    print('distance matrix is {} x {}; found {} matching hashes.'.format(mat.shape[0], mat.shape[1], len(hashes_to_tax)))

    print('clustering by togetherness & assigning taxonomy!')
    rootnode, nodelist, node_id_to_tax = do_cluster(mat, hashes_to_tax)

    if args.pickle_tree:
        print('pickling tree & taxonomy to file', args.pickle_tree)
        with open(args.pickle_tree, 'wb') as fp:
            dump((rootnode, nodelist, node_id_to_tax), fp)


if __name__ == '__main__':
    main()
