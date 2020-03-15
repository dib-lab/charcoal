#! /usr/bin/env python
"""
Take output of 'process-genome.py', match hashes against many metagenomes.
"""
import sys
import argparse
import sourmash
import screed
from pickle import load
import numpy as np
import csv
import os


def main():
    p = argparse.ArgumentParser()
    p.add_argument('hashes_pickle')
    p.add_argument('metagenome_sigs_list')
    p.add_argument('matrix_csv_out')
    p.add_argument('-d', '--metagenome-sigs-dir', default=None)
    args = p.parse_args()

    with open(args.hashes_pickle, 'rb') as fp:
        hash_to_lengths = load(fp)
    print('loaded {} hashes from {}'.format(len(hash_to_lengths), args.hashes_pickle))

    with open(args.metagenome_sigs_list, 'rt') as fp:
        metagenome_prefixes = [ x.strip() for x in fp ]

    if args.metagenome_sigs_dir:
        metagenome_prefixes = [ os.path.join(args.metagenome_sigs_dir, k) for k in metagenome_prefixes ]

    matrix = np.zeros((len(metagenome_prefixes), len(hash_to_lengths)), int)

    for i, prefix in enumerate(metagenome_prefixes):
        print('loading', prefix)
        ss = sourmash.load_one_signature(prefix, ksize=31)

        mins = ss.minhash.get_mins(with_abundance=True)

        m = 0
        for j, query in enumerate(sorted(hash_to_lengths)):
            count = mins.get(query, 0)
            if count:
                m += 1

                matrix[i][j] = count

        print('...', i, len(metagenome_prefixes), prefix, len(hash_to_lengths), m)

    print(len(metagenome_prefixes), len(hash_to_lengths))
    with open(args.matrix_csv_out, 'wt') as outfp:
        w = csv.writer(outfp)
        for i in range(len(metagenome_prefixes)):
            y = []
            for j in range(len(hash_to_lengths)):
                y.append('{}'.format(matrix[i][j]))
            w.writerow(y)

    return 0


if __name__ == '__main__':
    sys.exit(main())