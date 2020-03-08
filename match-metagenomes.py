#! /usr/bin/env python
import sys
import argparse
import sourmash
import screed
from pickle import load


def main():
    p = argparse.ArgumentParser()
    p.add_argument('hashes_pickle')
    p.add_argument('metagenome_sigs', nargs='+')
    args = p.parse_args()

    with open(args.hashes_pickle, 'rb') as fp:
        hash_to_lengths = load(fp)
    print('loaded {} hashes from {}'.format(len(hash_to_lengths), args.hashes_pickle))

    for metagenome in args.metagenome_sigs:
        ss = sourmash.load_one_signature(metagenome)

        mins = ss.minhash.get_mins()
        minset = set(mins)

        m = 0
        for query in hash_to_lengths:
            if query in minset:
                m += 1

        print('...', metagenome, len(hash_to_lengths), m)
            

    return 0


if __name__ == '__main__':
    sys.exit(main())
