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

import utils


def do_cluster(mat, hashes_to_tax):
    """
    Use scipy.cluster.hierarchy to cluster the matrix.
    """
    n_hashes = mat.shape[1]
    assert len(hashes_to_tax) == n_hashes

    # do the clustering...
    Y = sch.linkage(mat, method='complete')
    rootnode, nodelist = sch.to_tree(Y, rd=True)

    # now, track the taxonomy <-> hashval <-> cluster node.
    hashlist = list(sorted(hashes_to_tax))
    node_id_to_tax = {}

    # here we are assuming that the nodelist is in the same order as the
    # original labels. Is that right? Yes, it is. Just make sure the
    # hashes are sorted (same as when we created the matrix).
    for pos in range(n_hashes):
        node = nodelist[pos]
        assert node.get_id() == pos
        hashval = hashlist[pos]
        tax_info = hashes_to_tax[hashval]

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


def main():
    p = argparse.ArgumentParser()
    p.add_argument('matrix', help='output of match_metagenomes')
    p.add_argument('taxhashes', help='output of genome_shred_to_tax')
    p.add_argument('--pickle-tree', default=None)
    args = p.parse_args()

    # output of match_metagenomes
    print('calculating distance matrix from', args.matrix)
    with open(args.matrix, 'rb') as fp:
        matrix_obj = load(fp)
    mat, n_orig_hashes = utils.make_distance_matrix(matrix_obj.mat)

    # output of genome_shred_to_tax
    with open(args.taxhashes, 'rb') as fp:
        hashes_to_tax = load(fp)

    # some basic validation
    assert matrix_obj.ksize == hashes_to_tax.ksize
    print(matrix_obj.genome_file, hashes_to_tax.genome_file)
    assert matrix_obj.genome_file == hashes_to_tax.genome_file
    assert matrix_obj.query_hashlist == list(sorted(hashes_to_tax))
    assert matrix_obj.query_fragment_size == hashes_to_tax.fragment_size
    assert n_orig_hashes == len(hashes_to_tax), "mismatch! was same --scaled used to compute these?"
    assert mat.shape[1] == len(hashes_to_tax)

    print('distance matrix is {} x {}; found {} matching hashes.'.format(mat.shape[0], mat.shape[1], len(hashes_to_tax)))

    print('clustering by togetherness & assigning taxonomy!')
    rootnode, nodelist, node_id_to_tax = do_cluster(mat, hashes_to_tax)

    if args.pickle_tree:
        print('pickling tree & taxonomy to file', args.pickle_tree)
        with open(args.pickle_tree, 'wb') as fp:
            dump((rootnode, nodelist, node_id_to_tax), fp)


if __name__ == '__main__':
    main()
