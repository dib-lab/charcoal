#! /usr/bin/env python
"""
Assign taxonomy to shredded fragments in genomes.
"""
import sys
import argparse
from pickle import dump
import csv
from collections import defaultdict

import sourmash
import screed
from sourmash.lca import lca_utils
import utils                              # charcoal utils


def summarize(hashvals, dblist, threshold):
    """
    Classify 'hashvals' using the given list of databases.

    Insist on at least 'threshold' counts of a given lineage before taking
    it seriously.

    Return (lineage, counts) where 'lineage' is a tuple of LineagePairs.
    """

    # gather assignments from across all the databases
    assignments = lca_utils.gather_assignments(hashvals, dblist)

    # now convert to trees -> do LCA & counts
    counts = lca_utils.count_lca_for_assignments(assignments)

    # ok, we now have the LCAs for each hashval, and their number
    # of counts. Now aggregate counts across the tree, going up from
    # the leaves.
    aggregated_counts = defaultdict(int)
    for lca, count in counts.most_common():
        if count < threshold:
            break

        if not lca:
            aggregated_counts[lca] += count

        # climb from the lca to the root.
        aggregated_counts[lca] += count

    return aggregated_counts


def classify_signature(mh, db_list, threshold):
    # gather assignments from across all the databases
    assignments = lca_utils.gather_assignments(mh.get_mins(),
                                               db_list)

    # now convert to trees -> do LCA & counts
    counts = lca_utils.count_lca_for_assignments(assignments)

    # ok, we now have the LCAs for each hashval, and their number of
    # counts. Now build a tree across "significant" LCAs - those above
    # threshold.

    tree = {}

    for lca, count in counts.most_common():
        if count < threshold:
            break

        # update tree with this set of assignments
        lca_utils.build_tree([lca], tree)

    status = 'nomatch'
    if not tree:
        return [], status

    # now find lowest-common-ancestor of the resulting tree.
    lca, reason = lca_utils.find_lca(tree)
    if reason == 0:               # leaf node
        status = 'found'
    else:                         # internal node => disagreement
        status = 'disagree'

    return lca, status


def main():
    p = argparse.ArgumentParser()
    p.add_argument('lca_db')
    p.add_argument('genome')
    p.add_argument('output')
    p.add_argument('--fragment', default=100000, type=int)
    p.add_argument('--save-tax-hashes', default=None)
    args = p.parse_args()

    db, ksize, scaled = lca_utils.load_single_database(args.lca_db)
    mh_factory = sourmash.MinHash(n=0, ksize=ksize, scaled=scaled)
    print('**', ksize, scaled)

    n = 0
    m = 0
    sum_bp = 0
    sum_missed_bp = 0

    outfp = open(args.output, 'wt')
    w = csv.writer(outfp)
    w.writerow(['filename', 'contig', 'begin', 'end', 'lca', 'lca_rank', 'classified_as', 'classify_reason'])

    hashes_to_tax = utils.HashesToTaxonomy(args.genome,
                                           ksize,
                                           scaled,
                                           args.fragment,
                                           args.lca_db)

    #
    # iterate over all contigs in genome file
    #
    for record in screed.open(args.genome):
        # fragment longer contigs into smaller regions
        for start in range(0, len(record.sequence), args.fragment):
            seq = record.sequence[start:start + args.fragment]
            n += 1
            sum_bp += len(seq)

            # for each fragment, construct hashes
            mh = mh_factory.copy_and_clear()
            mh.add_sequence(seq, force=True)
            if not mh:
                sum_missed_bp += len(seq)
                continue

            # summarize & classify hashes; probably redundant code here.
            lineage_counts = summarize(mh.get_mins(), [db], 1)
            classify_lca, reason = classify_signature(mh, [db], 1)

            # output a CSV containing all of the lineage counts
            # (do we use this for anything?)
            for k in lineage_counts:
                lca_str = lca_utils.display_lineage(k, truncate_empty=False)
                classify_lca_str = lca_utils.display_lineage(classify_lca, truncate_empty=False)
                rank = ""
                if k:
                    rank = k[-1].rank
                w.writerow((args.genome, record.name, start, start + args.fragment, lca_str, rank, classify_lca_str, reason))

            # construct the hashes_to_tax dictionary from the minimum
            # of the hashes in the contig; this will match the
            # results from process_genome.
            min_of_mh = min(mh.get_mins())
            hashes_to_tax[min_of_mh] = classify_lca

            m += 1
            min_value = min(mh.get_mins())

    print('{} contigs / {} bp, {} hash values (missing {} contigs / {} bp)'.format(n, sum_bp, len(hashes_to_tax), n - m, sum_missed_bp))

    if args.save_tax_hashes:
        with open(args.save_tax_hashes, 'wb') as fp:
            dump(hashes_to_tax, fp)

    return 0


if __name__ == '__main__':
    sys.exit(main())
