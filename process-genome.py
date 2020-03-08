#! /usr/bin/env python
import sys
import argparse
import sourmash
import screed
from pickle import dump


def main():
    p = argparse.ArgumentParser()
    p.add_argument('genome')
    p.add_argument('output')
    p.add_argument('-k', '--ksize', default=31, type=int)
    p.add_argument('--scaled', default=1000, type=int)
    args = p.parse_args()

    mh_factory = sourmash.MinHash(n=0, ksize=args.ksize, scaled=args.scaled)

    hash_to_lengths = {}

    n = 0
    m = 0
    sum_bp = 0
    sum_missed_bp = 0
    for record in screed.open(args.genome):
        n += 1
        sum_bp += len(record.sequence)

        mh = mh_factory.copy_and_clear()
        mh.add_sequence(record.sequence, force=True)
        if not mh:
            sum_missed_bp += len(record.sequence)
            continue

        m += 1
        min_value = min(mh.get_mins())
        hash_to_lengths[min_value] = len(record.sequence)

    print('{} contigs / {} bp, {} hash values (missing {} contigs / {} bp)'.format(n, sum_bp, len(hash_to_lengths), n - m, sum_missed_bp))

    with open(args.output, 'wb') as fp:
        dump(hash_to_lengths, fp)
    
    return 0


if __name__ == '__main__':
    sys.exit(main())
