#! /usr/bin/env python
"""
Do gather matches on contigs => taxonomy, and save to JSON
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
from .utils import (gather_at_rank, get_ident, ContigGatherInfo)


def main(args):
    "Main entry point for scripting. Use cmdline for command line entry."
    genomebase = os.path.basename(args.genome)
    match_rank = args.match_rank

    # load taxonomy CSV
    tax_assign, _ = load_taxonomy_assignments(args.lineages_csv,
                                              start_column=2)
    print(f'loaded {len(tax_assign)} tax assignments.')

    # load the genome signature
    genome_sig = sourmash.load_one_signature(args.genome_sig)

    # load the matches from prefetch as a picklist
    picklist = sourmash.picklist.SignaturePicklist('prefetch',
                                                   pickfile=args.matches_csv)
    picklist.load(allow_empty=True)

    # load all of the matches in the database, as found by prefetch;
    # select on them; and then aggregate into MultiIndex.
    # CTB note: currently, this loads all the signatures into memory.
    # Alternatively we could do something with LazyLoadedIndex maybe?

    siglist = []
    for filename in args.databases:
        db = sourmash.load_file_as_index(filename)
        db = db.select(picklist=picklist)
        siglist += list(db.signatures())

    print(f"loaded {len(siglist)} matches from '{args.matches_csv}'")

    # Hack for examining members of our search database: remove exact matches.
    new_siglist = []
    for ss in siglist:
        if genome_sig.similarity(ss) == 1.0:
            print(f'removing an identical match: {ss.name}')
        else:
            new_siglist.append(ss)
    siglist = new_siglist

    # if, after removing exact match(es), there is nothing left, quit.
    # (but write an empty JSON file so that snakemake workflows don't
    # complain.)
    if not siglist:
        print('no non-identical matches for this genome, exiting.')
        contigs_tax = {}
        with open(args.json_out, 'wt') as fp:
            fp.write(json.dumps(contigs_tax))
        return 0

    # remove duplicate signatures in matches
    # workaround for issue of duplicate sigs in SBT, see sourmash/#1171
    new_siglist = []
    seen_md5 = set()
    for ss in siglist:
        ss_md5 = ss.md5sum()
        if not ss_md5 in seen_md5:
            new_siglist.append(ss)
            seen_md5.add(ss_md5)    
        else:
            print(f'removing a duplicate match: {ss.name}')
    siglist = new_siglist

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
        # look at each contig individually
        mh = empty_mh.copy_and_clear()
        mh.add_sequence(record.sequence, force=True)
        # collect all the gather results at genus level, together w/counts;
        # here, results is a list of (lineage, count) tuples.
        results = list(gather_at_rank(mh, lca_db, lin_db, match_rank))
        # store together with size of sequence.
        info = ContigGatherInfo(len(record.sequence), len(mh), results)
        contigs_tax[record.name] = info

    print(f"Processed {len(contigs_tax)} contigs.")

    # save!
    with open(args.json_out, 'wt') as fp:
        fp.write(json.dumps(contigs_tax))

    return 0


def cmdline(sys_args):
    "Command line entry point w/argparse action."
    p = argparse.ArgumentParser(sys_args)
    p.add_argument('--genome', help='genome file', required=True)
    p.add_argument('--genome-sig', help='genome sig', required=True)
    p.add_argument('--matches-csv', help='all relevant matches', required=True)
    p.add_argument('--databases', help='sourmash databases', required=True,
                   nargs='+')
    p.add_argument('--lineages-csv', help='lineage spreadsheet', required=True)
    p.add_argument('--force', help='continue past survivable errors',
                   action='store_true')

    p.add_argument('--json-out',
                   help='JSON-format output file of all tax results',
                   required=True)
    p.add_argument('--match-rank', required=True)
    args = p.parse_args()

    return main(args)


# execute this, when run with `python -m`.
if __name__ == '__main__':
    returncode = cmdline(sys.argv[1:])
    sys.exit(returncode)
