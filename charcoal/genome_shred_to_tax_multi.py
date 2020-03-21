#! /usr/bin/env python
"""
Assign taxonomy to shredded fragments in many genomes.

This does the same thing as genome_shred_to_tax, but for many genomes at once.
"""
import sys
import argparse
from pickle import dump
import csv
from collections import defaultdict
import os

import sourmash
import screed
from sourmash.lca import lca_utils
import utils                              # charcoal utils
from genome_shred_to_tax import summarize, classify_signature


def main():
    p = argparse.ArgumentParser()
    p.add_argument('lca_db')
    p.add_argument('genomes', nargs='+')
    p.add_argument('--csv-output-template', default=None)
    p.add_argument('--fragment', default=100000, type=int)
    p.add_argument('--save-tax-hashes-template', default=None)
    args = p.parse_args()

    assert args.csv_output_template

    db, ksize, scaled = lca_utils.load_single_database(args.lca_db)
    mh_factory = sourmash.MinHash(n=0, ksize=ksize, scaled=scaled)
    print('**', ksize, scaled)

    for genome in args.genomes:
        genome_base = os.path.basename(genome)
        output = args.csv_output_template.format(genome=genome_base)

        save_tax_hashes = None
        if args.save_tax_hashes_template:
            save_tax_hashes = args.save_tax_hashes_template.format(genome=genome_base)

        
        n = 0
        m = 0
        sum_bp = 0
        sum_missed_bp = 0

        outfp = open(output, 'wt')
        w = csv.writer(outfp)
        w.writerow(['filename', 'contig', 'begin', 'end', 'lca', 'lca_rank', 'classified_as', 'classify_reason'])

        hashes_to_tax = utils.HashesToTaxonomy(genome,
                                               ksize,
                                               scaled,
                                               args.fragment,
                                               args.lca_db)

        #
        # iterate over all contigs in genome file, fragmenting them.
        #
        shredder = utils.GenomeShredder(genome, args.fragment)
        for name, seq, start, end in shredder:
            n += 1
            sum_bp += len(seq)

            # for each fragment, construct hashes
            mh = mh_factory.copy_and_clear()
            mh.add_sequence(seq, force=True)
            if not mh:
                sum_missed_bp += len(seq)
                continue

            # summarize & classify hashes; probably redundant code here...
            lineage_counts = summarize(mh.get_mins(), [db], 1)
            classify_lca, reason = classify_signature(mh, [db], 1)

            # output a CSV containing all of the lineage counts
            # (do we use this for anything?)
            for k in lineage_counts:
                lca_str = lca_utils.display_lineage(k, truncate_empty=False)
                classify_lca_str = lca_utils.display_lineage(classify_lca,
                                                             truncate_empty=False)
                rank = ""
                if k:
                    rank = k[-1].rank
                w.writerow((genome, name, start, end,
                            lca_str, rank, classify_lca_str, reason))

            # construct the hashes_to_tax dictionary from the minimum
            # of the hashes in the contig; this will match the
            # results from process_genome.
            min_of_mh = min(mh.get_mins())
            hashes_to_tax[min_of_mh] = classify_lca

            m += 1
            min_value = min(mh.get_mins())

        # done! summarize to output.
        print('{} contigs / {} bp, {} hash values (missing {} contigs / {} bp)'.format(n, sum_bp, len(hashes_to_tax), n - m, sum_missed_bp))

        if save_tax_hashes:
            with open(save_tax_hashes, 'wb') as fp:
                dump(hashes_to_tax, fp)

    return 0


if __name__ == '__main__':
    sys.exit(main())
