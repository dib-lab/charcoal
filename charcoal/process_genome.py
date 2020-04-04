#! /usr/bin/env python
"""
Assign hashes to contigs and/or shredded fragments in genomes.
"""
import sys
import argparse
import sourmash
import screed
from pickle import dump

from . import utils                              # charcoal utils


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--genome', required=True)
    p.add_argument('--save-hashes', required=True)
    p.add_argument('-k', '--ksize', default=31, type=int)
    p.add_argument('--scaled', default=1000, type=int)
    p.add_argument('--fragment', default=0, type=int)
    p.add_argument('--stats', default=None)
    args = p.parse_args()

    mh_factory = sourmash.MinHash(n=0, ksize=args.ksize, scaled=args.scaled)

    hash_to_lengths = utils.HashesToLengths(args.genome,
                                            args.ksize, args.scaled,
                                            args.fragment)

    n = 0
    m = 0
    sum_bp = 0
    sum_missed_bp = 0

    statsfp = None
    if args.stats:
        statsfp = open(args.stats, 'wt')

    #
    # iterate over all contigs in genome file
    #
    shredder = utils.GenomeShredder(args.genome, args.fragment)
    for name, seq, start, end in shredder:
        n += 1
        sum_bp += len(seq)

        # for each contig, calculate scaled hashes
        mh = mh_factory.copy_and_clear()
        mh.add_sequence(seq, force=True)
        if statsfp and len(seq) == args.fragment:
            print('{}'.format(len(mh)), file=statsfp)

        # none assigned? so sad. record and move on.
        if not mh:
            sum_missed_bp += len(seq)
            continue

        # track the minimum of these for further analysis.
        m += 1
        min_value = min(mh.get_mins())
        hash_to_lengths[min_value] = len(seq)

    # some summary output
    print('{} contigs / {} bp, {} hash values (missing {} contigs / {} bp)'.format(n, sum_bp, len(hash_to_lengths), n - m, sum_missed_bp))

    # save pickled version
    with open(args.save_hashes, 'wb') as fp:
        dump(hash_to_lengths, fp)
    
    return 0


if __name__ == '__main__':
    sys.exit(main())
