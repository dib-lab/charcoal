#! /usr/bin/env python
"""
Output a more usable JSON format for the together tree.

Guide to mappings output in the JSON file:

    # hashval -> (contig name, fragment start, fragment end)
    output['hashes_to_fragment'] = hashes_to_fragment

    # leaf_id -> hashval
    output['leaves_to_hashval'] = leaves_to_hashval

    # hashval -> tax
    output['hashes_to_tax'] = hashes_to_tax

    # node_id -> tax
    output['node_id_to_tax'] = node_id_to_tax2

    # nodelist / cluster hierarchy
    output['node_to_children'] = node_to_children

"""
import sys
import argparse
import sourmash
from sourmash.lca import lca_utils
import screed
import json
from pickle import load
import csv


def main():
    p = argparse.ArgumentParser()
    p.add_argument('genome')
    p.add_argument('taxhashes', help='output of genome_shred_to_tax')
    p.add_argument('tree', help='output of combine_tax_togetherness')
    p.add_argument('json_output')
    p.add_argument('-k', '--ksize', default=31, type=int)
    p.add_argument('--scaled', default=1000, type=int)
    p.add_argument('--fragment', default=0, type=int)
    args = p.parse_args()

    assert args.fragment, "must specify --fragment"

    with open(args.tree, 'rb') as fp:
        (rootnode, nodelist, node_id_to_tax) = load(fp)

    # output of genome_shred_to_tax
    with open(args.taxhashes, 'rb') as fp:
        hashes_to_tax = load(fp)

    mh_factory = sourmash.MinHash(n=0, ksize=args.ksize, scaled=args.scaled)

    hashes_to_fragment = {}

    #
    # iterate over all contigs in genome file
    #
    for record in screed.open(args.genome):
        # fragment longer contigs into smaller regions?
        if args.fragment:
            for start in range(0, len(record.sequence), args.fragment):
                # fragment contigs
                seq = record.sequence[start:start + args.fragment]

                # for each contig, calculate scaled hashes
                mh = mh_factory.copy_and_clear()
                mh.add_sequence(seq, force=True)

                minset = set(mh.get_mins())

                if minset:
                    min_hash_val = min(minset)

                    hashes_to_fragment[min_hash_val] = (record.name,
                                                        start,
                                                        start + len(seq))

        # deal directly with contigs in the genome file
        else:
            assert 0
            n += 1

            mh = mh_factory.copy_and_clear()
            mh.add_sequence(record.sequence, force=True)
            if not mh:
                sum_missed_bp += len(record.sequence)
                continue

            m += 1
            min_value = min(mh.get_mins())
            hash_to_lengths[min_value] = len(record.sequence)

    # construct mapping from leaf ID to hashval
    leaves_to_hashval = {}
    hashes = sorted(hashes_to_fragment)
    for leaf_id in range(len(hashes_to_fragment)):
        leaves_to_hashval[leaf_id] = hashes[leaf_id]

    # nodes to children IDs
    node_to_children = {}
    for pos, node in enumerate(nodelist):
        node_id = node.get_id()
        l_id = -1
        r_id = -1
        if not node.is_leaf():
            l_id = node.get_left().get_id()
            r_id = node.get_right().get_id()

        node_to_children[node_id] = (l_id, r_id)

    # node ID to tax, serializable
    node_id_to_tax2 = {}
    for k, v in node_id_to_tax.items():
        node_id_to_tax2[k] = list([ list(x) for x in v ])

    #
    # ok, our json file will contain the following:
    #
    output = {}
    
    # hashval -> (contig name, fragment start, fragment end)
    output['hashes_to_fragment'] = hashes_to_fragment

    # leaf_id -> hashval
    output['leaves_to_hashval'] = leaves_to_hashval

    # hashval -> tax
    output['hashes_to_tax'] = hashes_to_tax.d

    # node_id -> tax
    output['node_id_to_tax'] = node_id_to_tax2

    # nodelist / cluster hierarchy
    output['node_to_children'] = node_to_children

    with open(args.json_output, 'wt') as outfp:
        outfp.write(json.dumps(output))

    with open(args.json_output + '.hashes_to_fragment.csv', 'wt') as outfp:
        w = csv.writer(outfp)

        w.writerow(['hashval', 'contig', 'start', 'end'])
        for hashval, (contig, start, end) in hashes_to_fragment.items():
            w.writerow([hashval, contig, start, end])

    with open(args.json_output + '.leaves_to_hashval.csv', 'wt') as outfp:
        w = csv.writer(outfp)

        w.writerow(['leaf_id', 'hashval'])
        for leaf_id, hashval in leaves_to_hashval.items():
            w.writerow([leaf_id, hashval])

    with open(args.json_output + '.hashes_to_tax.csv', 'wt') as outfp:
        w = csv.writer(outfp)

        w.writerow(['hashval', 'lineage'])
        for hashval, lca in hashes_to_tax.items():
            lineage_str = lca_utils.display_lineage(lca, truncate_empty=False, include_strain=True)
            x = [hashval, lineage_str]
            w.writerow(x)

    with open(args.json_output + '.node_id_to_tax.csv', 'wt') as outfp:
        w = csv.writer(outfp)

        w.writerow(['node_id'])
        for node_id, taxset in node_id_to_tax2.items():
            x = [node_id]
            for tax in taxset:
                x.append(lca_utils.display_lineage(tax, truncate_empty=False, include_strain=True))
            w.writerow(x)

    with open(args.json_output + '.node_to_children.csv', 'wt') as outfp:
        w = csv.writer(outfp)

        w.writerow(['node_id', 'left_child_id', 'right_child_id'])
        for node_id, (lid, rid) in node_to_children.items():
            w.writerow([node_id, lid, rid])

    return 0


if __name__ == '__main__':
    sys.exit(main())
