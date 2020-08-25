#! /usr/bin/env python
"""
Remove bad contigs based solely on taxonomy.
"""
import sys
import argparse
import gzip
from collections import Counter, defaultdict
import csv
import os.path
import json
import glob

import screed

import sourmash
from sourmash.lca.command_index import load_taxonomy_assignments
from sourmash.lca import LCA_Database, LineagePair

from . import utils
from . import lineage_db
from .lineage_db import LineageDB
from .version import version
from .utils import (get_idents_for_hashval, gather_lca_assignments,
    count_lca_for_assignments, pretty_print_lineage, pretty_print_lineage2,
    WriteAndTrackFasta, gather_at_rank, get_ident, summarize_at_rank)


GATHER_MIN_MATCHES=3
F_IDENT_THRESHOLD=0.1
F_MAJOR_THRESHOLD=0.2


from enum import Enum
class ContigInfo(Enum):
    "enum for contig information"
    CLEAN = 1                             # proven clean
    DIRTY = 2                             # proven dirty
    NO_IDENT = 3                          # no identified hashes in contig
    NO_HASH = 4                           # no hashes in this contig


def kb(bp):
    return int(bp/1000)



def calculate_contam(info_obj, rank, filter_names=None):
    genome_lin = info_obj.genome_tax
    good_names = set()
    good_n = 0
    good_bp = 0
    bad_names = set()
    bad_n = 0
    bad_bp = 0

    for contig_name, (contig_len, contig_taxlist) in info_obj.contigs_d.items():
        if filter_names and contig_name in filter_names:
            continue

        contig_taxlist_at_rank = summarize_at_rank(contig_taxlist, rank)
        top_hit = None
        if contig_taxlist_at_rank:
            top_hit = contig_taxlist_at_rank[0][0]

        if genome_lin and top_hit and not utils.is_lineage_match(genome_lin, top_hit, rank):
            bad_names.add(contig_name)
            bad_n += 1
            bad_bp += contig_len
        else:
            good_names.add(contig_name)
            good_n += 1
            good_bp += contig_len

    return (good_names, good_n, good_bp, bad_names, bad_n, bad_bp)


class GenomeAndContigsInfo(object):
    def __init__(self, genome_name, genome_tax, contigs_d):
        self.genome_name = genome_name
        self.genome_tax = genome_tax
        self.contigs_d = contigs_d


class ContigReport(object):
    # output in addition: genomefile, genome taxonomy, match_rank
    def __init__(self, contig_name, bp, info, reason, lineage, total_hashes,
                 ident_hashes, match_hashes):
        self.contig_name = contig_name
        self.bp = bp
        self.info = info
        self.reason = reason
        self.lineage = lineage
        self.total_hashes = total_hashes
        self.ident_hashes = ident_hashes
        self.match_hashes = match_hashes


def guess_tax_by_gather(entire_mh, lca_db, lin_db, match_rank, report_fp):
    "Guess likely taxonomy using gather."
    sum_ident = 0
    first_lin = ()
    first_count = 0

    for lin, count in gather_at_rank(entire_mh, lca_db, lin_db, match_rank):
        if count >= GATHER_MIN_MATCHES:
            # record the first lineage we come across as likely lineage.
            if not first_lin:
                first_lin = lin
                first_count = count

        sum_ident += count

    f_ident = sum_ident / len(entire_mh)
    f_major = first_count / sum_ident

    return first_lin, f_ident, f_major


def choose_genome_lineage(guessed_genome_lineage, provided_lineage, match_rank,
                          f_ident, f_major, report):

    comment = ""
    genome_lineage = None

    if provided_lineage:
        if utils.is_lineage_match(provided_lineage, guessed_genome_lineage, match_rank):
            report(f'(provided lineage agrees with k-mer classification at {match_rank} level)')
        elif guessed_genome_lineage:
            report(f'(provided lineage disagrees with k-mer classification at or above {match_rank} level)')
        else:
            pass

        genome_lineage = utils.pop_to_rank(provided_lineage, match_rank)
        report(f'\nUsing provided lineage as genome lineage.')
    else:
        if f_ident < F_IDENT_THRESHOLD:
            report(f'** ERROR: fraction of total identified hashes (f_ident) < {F_IDENT_THRESHOLD*100:.0f}%.')
            comment = f"too few identifiable hashes; f_ident < {F_IDENT_THRESHOLD*100:.0f}%. provide a lineage for this genome."
        elif f_major < F_MAJOR_THRESHOLD:
            report(f'** ERROR: fraction of identified hashes in major lineage (f_major) < {F_MAJOR_THRESHOLD*100:.0f}%.')
            comment = f"too few hashes in major lineage; f_major < {F_MAJOR_THRESHOLD*100:.0f}%. provide a lineage for this genome."
        else:
            genome_lineage = utils.pop_to_rank(guessed_genome_lineage, match_rank)
            report(f'Using majority gather lineage as genome lineage.')

    return genome_lineage, comment


def get_genome_taxonomy(matches_filename, genome_sig_filename, provided_lineage,
                        tax_assign, match_rank):
    with open(matches_filename, 'rt') as fp:
        try:
            siglist = list(sourmash.load_signatures(fp, do_raise=True, quiet=True))
        except sourmash.exceptions.SourmashError:
            siglist = None

    if not siglist:
        print('no matches for this genome, exiting.')
        return None, 'no matches for this genome, exiting.', 0.0, 0.0

    # construct a template minhash object that we can use to create new 'uns
    empty_mh = siglist[0].minhash.copy_and_clear()
    ksize = empty_mh.ksize
    scaled = empty_mh.scaled
    moltype = empty_mh.moltype

    genome_sig = sourmash.load_one_signature(genome_sig_filename)
    entire_mh = genome_sig.minhash

    assert entire_mh.scaled == scaled

    # Hack for examining members of our search database: remove exact matches.
    new_siglist = []
    identical_match_removed = False
    for ss in siglist:
        if entire_mh.similarity(ss.minhash) < 1.0:
            new_siglist.append(ss)
        else:
            if provided_lineage and provided_lineage != 'NA':
                print(f'found exact match: {ss.name()}. removing.')
                identical_match_removed = True
            else:
                print(f'found exact match: {ss.name()}. but no provided lineage! exiting.')
                return None, f'found exact match: {ss.name()}. but no provided lineage! exiting.', 1.0, 1.0

    # ...but leave exact matches in if they're the only matches, I guess!
    if new_siglist:
        siglist = new_siglist

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

    # calculate lineage from majority vote on LCA
    guessed_genome_lineage, f_major, f_ident = \
         guess_tax_by_gather(entire_mh, lca_db, lin_db, match_rank, sys.stdout)

    print(f'Gather classification on this genome yields: {pretty_print_lineage(guessed_genome_lineage)}')

    if f_major == 1.0 and f_ident == 1.0:
        comment = "All genome hashes belong to one lineage! Nothing to do."
        return None, comment, f_major, f_ident

    # did we get a passed-in lineage assignment?
    provided_lin = ""
    if provided_lineage:
        provided_lin = [ LineagePair(rank, name) for (rank, name) in zip(sourmash.lca.taxlist(), provided_lin) if name.strip() ]
        print(f'Provided lineage from command line:\n   {sourmash.lca.display_lineage(provided_lin)}')

    # choose between the lineages
    genome_lineage, comment = choose_genome_lineage(guessed_genome_lineage,
                                                    provided_lin,
                                                    match_rank,
                                                    f_ident, f_major,
                                                    print)

    return genome_lineage, comment, f_major, f_ident


###

def main(args):
    "Main entry point for scripting. Use cmdline for command line entry."
    match_rank = args.match_rank
    tax_assign, _ = load_taxonomy_assignments(args.lineages_csv,
                                              start_column=3)
    print(f'loaded {len(tax_assign)} tax assignments.')

    genome_names = [ x.strip() for x in open(args.genome_list_file) ]
    genome_names = [ os.path.basename(n) for n in genome_names ]
    dirname = args.input_directory

    provided_lineages = {}
    if args.provided_lineages:
        with open(args.provided_lineages, 'rt') as fp:
            r = csv.reader(fp)
            for row in r:
                genome_name = row[0]
                lin = row[1:]
                provided_lineages[genome_name] = lin

    print(f"loaded {len(provided_lineages)} provided lineages")

    # process every genome individually.
    info_d = {}
    summary_d = {}
    fp = open(args.output, 'wt')
    summary_w = csv.writer(fp)
    summary_w.writerow(['genome', 'filter_at', 'override_filter_at', 'f_ident', 'f_major', 'total_bad_bp', 'superkingdom_bad_bp', 'phylum_bad_bp', 'class_bad_bp', 'order_bad_bp', 'family_bad_bp', 'genus_bad_bp', 'lineage', 'comment'])
    for genome_name in genome_names:
        matches_filename = os.path.join(dirname, genome_name + '.gather-matches.sig')
        genome_sig = os.path.join(dirname, genome_name + '.sig')
        lineage = provided_lineages.get(genome_name, '')
        contigs_json = os.path.join(dirname, genome_name + '.contigs-tax.json')

        genome_lineage, comment, f_major, f_ident = get_genome_taxonomy(matches_filename,
                                                      genome_sig,
                                                      lineage,
                                                      tax_assign, match_rank)
        if genome_lineage:
            filter_at = match_rank
        else:
            genome_lineage = []
            filter_at = 'none'

        # load contigs JSON file
        with open(contigs_json, 'rt') as fp:
            contigs_d = json.load(fp)
            for k in contigs_d:
                (size, v) = contigs_d[k]
                vv = []
                for (lin, count) in v:
                    vv.append(([ LineagePair(*x) for x in lin ], count))
                contigs_d[k] = (size, vv)

        # make an info object to store everything; probably not needed?
        info_obj = GenomeAndContigsInfo(genome_name, genome_lineage, contigs_d)

        info_d[genome_name] = info_obj

        eliminate = set()
        print(f'examining {genome_name} for contamination:')

        summary_d[genome_name] = {}

        total_bad_n = 0
        total_bad_bp = 0
        for rank in sourmash.lca.taxlist():
            x = calculate_contam(info_obj, rank, filter_names=eliminate)
            (good_names, good_n, good_bp, bad_names, bad_n, bad_bp) = x
            eliminate.update(bad_names)
            total_bad_n += bad_n
            total_bad_bp += bad_bp

            print(f'   {rank}: {bad_n} contigs w/ {kb(bad_bp)}kb')
            summary_d[genome_name][rank] = total_bad_bp
            if rank == match_rank:
                break

        print(f'   (total): {total_bad_n} contigs w/ {kb(total_bad_bp)}kb')
        summary_w.writerow([genome_name, filter_at, '', f'{f_ident:.03}', f'{f_major:.03}', total_bad_bp, summary_d[genome_name]['superkingdom'], summary_d[genome_name]['phylum'], summary_d[genome_name]['class'], summary_d[genome_name]['order'], summary_d[genome_name]['family'], summary_d[genome_name]['genus'], sourmash.lca.display_lineage(genome_lineage), comment])

    ####

    print(f"processed {len(info_d)} genomes.")

    return 0


def cmdline(sys_args):
    "Command line entry point w/argparse action."
    p = argparse.ArgumentParser(sys_args)
    p.add_argument('--input-directory', required=True)
    p.add_argument('--genome-list-file', required=True)
    p.add_argument('-o', '--output', required=True)
    #p.add_argument('--genome_sig', help='genome signature', required=True)
    p.add_argument('--lineages_csv', help='lineage spreadsheet', required=True)
    p.add_argument('--provided_lineages', help='provided lineages', required=True)
    #p.add_argument('--matches_sig', help='all relevant matches', required=True)
    #p.add_argument('--contigs_json', help='output of contigs_search', required=True)
    #p.add_argument('--summary', help='CSV one line output')
    #p.add_argument('--force', help='continue past survivable errors',
    #               action='store_true')

    #p.add_argument('--lineage', help=';-separated lineage down to genus level',
    #               default='NA')        # default is str NA
    p.add_argument('--match-rank', help='rank below which matches are _not_ contaminants', default='genus')
    #p.add_argument('--contig-report', help='contig report (CSV)')
    args = p.parse_args()

    main(args)


# execute this, when run with `python -m`.
if __name__ == '__main__':
    returncode = cmdline(sys.argv[1:])
    sys.exit(0)
