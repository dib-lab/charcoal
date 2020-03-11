#! /usr/bin/env python
import sys
import argparse
import sourmash
import screed
from pickle import dump
import csv
from sourmash.lca import lca_utils
import pprint
from collections import defaultdict


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


def main():
    p = argparse.ArgumentParser()
    p.add_argument('lca_db')
    p.add_argument('genome', nargs='+')
    p.add_argument('output')
    p.add_argument('--fragment', default=100000, type=int)
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
    w.writerow(['filename', 'contig', 'begin', 'end', 'lca', 'lca_rank'])

    #
    # iterate over all contigs in genome file
    #
    for genome in args.genome:
        for record in screed.open(genome):
            # fragment longer contigs into smaller regions?
            for start in range(0, len(record.sequence), args.fragment):
                seq = record.sequence[start:start + args.fragment]
                n += 1
                sum_bp += len(seq)

                mh = mh_factory.copy_and_clear()
                mh.add_sequence(seq, force=True)
                if not mh:
                    sum_missed_bp += len(seq)
                    continue

                lineage_counts = summarize(mh.get_mins(), [db], 1)

                for k in lineage_counts:
                    lca = lca_utils.display_lineage(k, truncate_empty=False)
                    try:
                        lca_rank = k[-1].rank
                    except IndexError:
                        lca_rank = "none"
                    w.writerow((genome, record.name, start, start + args.fragment, lca, lca_rank))

                m += 1
                min_value = min(mh.get_mins())

    return 0


if __name__ == '__main__':
    sys.exit(main())
