#! /usr/bin/env python
"""
Do gather matches on contigs.
"""
import sys
import argparse
import os.path
import json

import screed

import sourmash
from sourmash.lca.command_index import load_taxonomy_assignments
from sourmash.lca import LCA_Database

from .lineage_db import LineageDB
from .version import version
from .utils import (gather_at_rank, get_ident)


def main(args):
    "Main entry point for scripting. Use cmdline for command line entry."
    genomebase = os.path.basename(args.genome)
    match_rank = args.match_rank

    tax_assign, _ = load_taxonomy_assignments(args.lineages_csv,
                                              start_column=3)
    print(f'loaded {len(tax_assign)} tax assignments.')

    with open(args.matches_sig, 'rt') as fp:
        siglist = list(sourmash.load_signatures(fp))

    genome_sig = sourmash.load_one_signature(args.genome_sig)

    # Hack for examining members of our search database: remove exact matches.
    new_siglist = []
    for ss in siglist:
        if genome_sig.similarity(ss) < 1.0:
            new_siglist.append(ss)

    siglist = new_siglist

    if not siglist:
        print('no non-identical matches for this genome, exiting.')
        contigs_tax = {}
        with open(args.json_out, 'wt') as fp:
            fp.write(json.dumps(contigs_tax))
        return 0

    # construct a template minhash object that we can use to create new 'uns
    empty_mh = siglist[0].minhash.copy_and_clear()
    ksize = empty_mh.ksize
    scaled = empty_mh.scaled
    moltype = empty_mh.moltype

    # create empty LCA database to populate...
    lca_db = LCA_Database(ksize=ksize, scaled=scaled, moltype=moltype)
    lin_db = LineageDB()

    # ...with specific matches.
    for ss in siglist:
        ident = get_ident(ss)
        lineage = tax_assign[ident]

        lca_db.insert(ss, ident=ident)
        lin_db.insert(ident, lineage)

    print(f'loaded {len(siglist)} signatures & created LCA Database')

    print('')
    print(f'reading contigs from {genomebase}')

    screed_iter = screed.open(args.genome)
    contigs_tax = {}
    for n, record in enumerate(screed_iter):
        # look at each contig
        mh = empty_mh.copy_and_clear()
        mh.add_sequence(record.sequence, force=True)

        results = list(gather_at_rank(mh, lca_db, lin_db, match_rank))

        contigs_tax[record.name] = (len(record.sequence), results)

    print(f"Processed {len(contigs_tax)} contigs.")

    with open(args.json_out, 'wt') as fp:
        fp.write(json.dumps(contigs_tax))

    return 0


def cmdline(sys_args):
    "Command line entry point w/argparse action."
    p = argparse.ArgumentParser(sys_args)
    p.add_argument('--genome', help='genome file', required=True)
    p.add_argument('--genome_sig', help='genome sig', required=True)
    p.add_argument('--matches_sig', help='all relevant matches', required=True)
    p.add_argument('--lineages_csv', help='lineage spreadsheet', required=True)
    p.add_argument('--force', help='continue past survivable errors',
                   action='store_true')

    p.add_argument('--match-rank', help='rank below which matches are _not_ contaminants', default='genus')

    p.add_argument('--json-out', help='JSON-format output file of all tax results')
    args = p.parse_args()

    main(args)


# execute this, when run with `python -m`.
if __name__ == '__main__':
    returncode = cmdline(sys.argv[1:])
    sys.exit(0)
