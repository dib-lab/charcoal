#! /usr/bin/env python
"""
Do gather matches on contigs => taxonomy, and save to JSON
"""
import sys
import argparse
import os.path
import yaml
from collections import defaultdict

import screed

import sourmash
from sourmash.lca.command_index import load_taxonomy_assignments
from sourmash.lca import LCA_Database

from .lineage_db import LineageDB
from .version import version
from . import utils
from .utils import (get_ident, CSV_DictHelper, make_lineage)
from .compare_taxonomy import GATHER_MIN_MATCHES


def get_matches(mh, lca_db, lin_db, match_rank, threshold_bp):
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
        common = match_sig.minhash.count_common(query_sig.minhash)

        # retrieve identity
        match_ident = get_ident(match_sig)
        match_lineage = lin_db.ident_to_lineage[match_ident]
        
        yield match_ident, match_lineage, common
        
        # finish out gather algorithm!
        minhash.remove_many(match_sig.minhash.hashes)
        query_sig = sourmash.SourmashSignature(minhash)

    return accs


def main(args):
    "Main entry point for scripting. Use cmdline for command line entry."
    genomebase = os.path.basename(args.genome)

    # load hitlist
    hitlist = CSV_DictHelper(args.hitlist, 'genome')
    hitlist_entry = hitlist[genomebase]
    match_rank = hitlist_entry.filter_at
    if hitlist_entry.override_filter_at:
        match_rank = hitlist_entry.override_filter_at

    assert match_rank in ('superkingdom', 'phylum', 'class', 'order',
                          'family', 'genus'), match_rank

    # genome lineage:
    genome_lin = make_lineage(hitlist_entry.lineage)

    # load taxonomy CSV
    tax_assign, _ = load_taxonomy_assignments(args.lineages_csv,
                                              start_column=3)
    print(f'loaded {len(tax_assign)} tax assignments.')

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
        with open(args.yaml_out, 'wt') as fp:
            pass
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

    matches_info = {}
    matches_counts = defaultdict(int)

    screed_iter = screed.open(args.genome)
    threshold_bp = empty_mh.scaled * GATHER_MIN_MATCHES
    n = - 1
    for n, record in enumerate(screed_iter):
        # look at each contig individually
        mh = empty_mh.copy_and_clear()
        mh.add_sequence(record.sequence, force=True)

        # collect all the gather results at genus level, together w/counts;
        # here, results is a list of (lineage, count) tuples.
        results = list(utils.gather_at_rank(mh, lca_db, lin_db, match_rank))

        for acc, match_lin, count in get_matches(mh, lca_db, lin_db, match_rank,
                                          threshold_bp):
            # dirty match
            if not utils.is_lineage_match(genome_lin, match_lin, match_rank):
                if acc in matches_info:
                    assert matches_info[acc][0] == 'dirty'
                matches_info[acc] = ['dirty', utils.display_lineage(match_lin)]
            else:                     # clean
                if acc in matches_info:
                    assert matches_info[acc][0] == 'clean'
                matches_info[acc] = ['clean', utils.display_lineage(match_lin)]

            matches_counts[acc] += count
                    

    print(f"Processed {n + 1} contigs.")

    # save!
    with open(args.yaml_out, 'wt') as fp:
        out_dict = {}
        info_dict = {}
        info_dict['genome'] = genomebase
        info_dict['genome_lineage'] = utils.display_lineage(genome_lin)
        out_dict['query_info'] = info_dict

        matches_info_out = {}
        for acc, (match_type, lineage) in matches_info.items():
            acc_info = {}
            acc_info['lineage'] = lineage
            acc_info['match_type'] = match_type
            acc_info['counts'] = matches_counts[acc]
            matches_info_out[acc] = acc_info
        out_dict['matches'] = matches_info_out

        yaml.dump(out_dict, fp)

    return 0


def cmdline(sys_args):
    "Command line entry point w/argparse action."
    p = argparse.ArgumentParser(sys_args)
    p.add_argument('--genome', help='genome file', required=True)
    p.add_argument('--genome-sig', help='genome sig', required=True)
    p.add_argument('--matches-sig', help='all relevant matches', required=True)
    p.add_argument('--lineages-csv', help='lineage spreadsheet', required=True)
    p.add_argument('--hitlist', help='hitlist spreadsheet', required=True)
    p.add_argument('--force', help='continue past survivable errors',
                   action='store_true')

    p.add_argument('--yaml-out',
                   help='YAML-format output file of all matches',
                   required=True)
    args = p.parse_args()

    return main(args)


# execute this, when run with `python -m`.
if __name__ == '__main__':
    returncode = cmdline(sys.argv[1:])
    sys.exit(returncode)
