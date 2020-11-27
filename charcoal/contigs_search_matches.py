#! /usr/bin/env python
"""
Do gather matches on contigs => match names.
"""
import sys
import argparse
import os.path
import json

import screed

import sourmash
from sourmash.lca import LCA_Database

from .lineage_db import LineageDB
from .version import version
from .utils import (get_ident, CSV_DictHelper)
from .compare_taxonomy import GATHER_MIN_MATCHES

def get_matches(mh, lca_db, match_rank, threshold_bp):
    import copy
    minhash = copy.copy(mh)
    query_sig = sourmash.SourmashSignature(minhash)

    accs = set()

    # do the gather:
    while 1:
        results = lca_db.gather(query_sig, threshold_bp=threshold_bp)
        if not results:
            break

        (match, match_sig, _) = results[0]

        # retrieve identity
        match_ident = get_ident(match_sig)
        accs.add(match_ident)
        
        # finish out gather algorithm!
        minhash.remove_many(match_sig.minhash.hashes)
        query_sig = sourmash.SourmashSignature(minhash)

    return accs


def main(args):
    "Main entry point for scripting. Use cmdline for command line entry."
    genomebase = os.path.basename(args.genome)
    match_rank = 'genus'

    # load the genome signature
    genome_sig = sourmash.load_one_signature(args.genome_sig)

    # load all of the matches from search --containment in the database
    with open(args.matches_sig, 'rt') as fp:
        try:
            siglist = list(sourmash.load_signatures(fp, do_raise=True,
                                                    quiet=False))
        except sourmash.exceptions.SourmashError:
            siglist = []
    print(f"loaded {len(siglist)} matches from '{args.matches_sig}'")

    # Hack for examining members of our search database: remove exact matches.
    new_siglist = []
    for ss in siglist:
        if genome_sig.similarity(ss) == 1.0:
            print(f'removing an identical match: {ss.name}')
        else:
            new_siglist.append(ss)
    siglist = new_siglist

    # if, after removing exact match(es), there is nothing left, quit.
    # (but write an empty output file so that snakemake workflows don't
    # complain.)
    if not siglist:
        print('no non-identical matches for this genome, exiting.')
        with open(args.output, 'wt') as fp:
            fp.write('')
        return 0

    # construct a template minhash object that we can use to create new 'uns
    empty_mh = siglist[0].minhash.copy_and_clear()
    ksize = empty_mh.ksize
    scaled = empty_mh.scaled
    moltype = empty_mh.moltype

    # create empty LCA database to populate...
    lca_db = LCA_Database(ksize=ksize, scaled=scaled, moltype=moltype)

    # ...with specific matches.
    for ss in siglist:
        ident = get_ident(ss)
        lca_db.insert(ss, ident=ident)

    print(f'loaded {len(siglist)} signatures & created LCA Database')

    print('')
    print(f'reading contigs from {genomebase}')

    screed_iter = screed.open(args.genome)
    accs = set()
    threshold_bp = GATHER_MIN_MATCHES * empty_mh.scaled
    n = -1
    for n, record in enumerate(screed_iter):
        # look at each contig individually
        mh = empty_mh.copy_and_clear()
        mh.add_sequence(record.sequence, force=True)

        # collect all the accessions from gather
        accs.update(get_matches(mh, lca_db, match_rank, threshold_bp))

    print(f"Processed {n+1} contigs.")

    # save!
    with open(args.output, 'wt') as fp:
        fp.write("\n".join(accs) + "\n")

    return 0


def cmdline(sys_args):
    "Command line entry point w/argparse action."
    p = argparse.ArgumentParser(sys_args)
    p.add_argument('--genome', help='genome file', required=True)
    p.add_argument('--genome-sig', help='genome sig', required=True)
    p.add_argument('--matches-sig', help='all relevant matches', required=True)
    p.add_argument('--force', help='continue past survivable errors',
                   action='store_true')

    p.add_argument('--output',
                   help='list of match accessions',
                   required=True)
    args = p.parse_args()

    return main(args)


# execute this, when run with `python -m`.
if __name__ == '__main__':
    returncode = cmdline(sys.argv[1:])
    sys.exit(returncode)
