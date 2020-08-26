#! /usr/bin/env python
"""
Create a "hit list" of how much will be removed at what ranks.
"""
import sys
import argparse
import csv
import os.path

import sourmash
from sourmash.lca.command_index import load_taxonomy_assignments
from sourmash.lca import LCA_Database, LineagePair

from . import utils
from .lineage_db import LineageDB
from .utils import (gather_at_rank, get_ident, summarize_at_rank,
                    pretty_print_lineage, load_contigs_gather_json)


GATHER_MIN_MATCHES=3
F_IDENT_THRESHOLD=0.1
F_MAJOR_THRESHOLD=0.2


def kb(bp):
    return int(bp/1000)


def calculate_contam(genome_lin, contigs_d, rank, filter_names=None):
    good_names = dict()
    bad_names = dict()

    for contig_name, (contig_len, contig_taxlist) in contigs_d.items():
        if filter_names and contig_name in filter_names:
            continue

        contig_taxlist_at_rank = summarize_at_rank(contig_taxlist, rank)
        top_hit = None
        if contig_taxlist_at_rank:
            top_hit, count = contig_taxlist_at_rank[0]
            if count < GATHER_MIN_MATCHES:
                top_hit = None

        is_bad = False
        if genome_lin and top_hit and not utils.is_lineage_match(genome_lin, top_hit, rank):
            is_bad = True

            # rescue?
            for hit, count in contig_taxlist_at_rank[1:]:
                if utils.is_lineage_match(genome_lin, hit, rank):
                    is_bad = False

        if is_bad:
            bad_names[contig_name] = contig_len
        else:
            good_names[contig_name] = contig_len

    return (good_names, bad_names)


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
                          f_ident, f_major, min_f_ident, min_f_major, report):

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
        if f_ident < min_f_ident:
            report(f'** ERROR: fraction of total identified hashes (f_ident) < {min_f_ident*100:.0f}%.')
            comment = f"too few identifiable hashes; f_ident < {min_f_ident*100:.0f}%. provide a lineage for this genome."
        elif f_major < min_f_major:
            report(f'** ERROR: fraction of identified hashes in major lineage (f_major) < {min_f_major*100:.0f}%.')
            comment = f"too few hashes in major lineage; f_major < {min_f_major*100:.0f}%. provide a lineage for this genome."
        else:
            genome_lineage = utils.pop_to_rank(guessed_genome_lineage, match_rank)
            report(f'Using majority gather lineage as genome lineage.')

    return genome_lineage, comment


def get_genome_taxonomy(matches_filename, genome_sig_filename, provided_lineage,
                        tax_assign, match_rank, min_f_ident, min_f_major):
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
    for ss in siglist:
        if entire_mh.similarity(ss.minhash) < 1.0:
            new_siglist.append(ss)
        else:
            if provided_lineage and provided_lineage != 'NA':
                print(f'found exact match: {ss.name()}. removing.')
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
        provided_lin = [ LineagePair(rank, name) for (rank, name) in zip(sourmash.lca.taxlist(), provided_lineage) if name.strip() ]
        print(f'Provided lineage from command line:\n   {sourmash.lca.display_lineage(provided_lin)}')

    # choose between the lineages
    genome_lineage, comment = choose_genome_lineage(guessed_genome_lineage,
                                                    provided_lin,
                                                    match_rank,
                                                    f_ident, f_major,
                                                    min_f_ident, min_f_major,
                                                    print)

    return genome_lineage, comment, f_major, f_ident

###

def main(args):
    "Main entry point for scripting. Use cmdline for command line entry."
    match_rank = 'genus'

    # load taxonomy assignments for all the things
    tax_assign, _ = load_taxonomy_assignments(args.lineages_csv,
                                              start_column=3)
    print(f'loaded {len(tax_assign)} tax assignments.')

    # load in all the genome names
    genome_names = set([ x.strip() for x in open(args.genome_list_file) ])
    genome_names = set([ os.path.basename(n) for n in genome_names ])
    dirname = args.input_directory

    print(f"loaded list of {len(genome_names)} genomes.")

    # load the provided lineages file
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
    summary_d = {}
    for genome_name in genome_names:
        matches_filename = os.path.join(dirname, genome_name + '.matches.sig')
        genome_sig = os.path.join(dirname, genome_name + '.sig')
        lineage = provided_lineages.get(genome_name, '')
        contigs_json = os.path.join(dirname, genome_name + '.contigs-tax.json')

        x = get_genome_taxonomy(matches_filename,
                                genome_sig,
                                lineage,
                                tax_assign, match_rank,
                                args.min_f_ident,
                                args.min_f_major)
        genome_lineage, comment, f_major, f_ident = x

        # did we get a lineage for this genome? if so, propose filtering at
        # default rank 'match_rank', otherwise ...do not filter.
        if genome_lineage:
            filter_at = match_rank
        else:
            genome_lineage = []
            filter_at = 'none'

        # load contigs tax
        contigs_d = load_contigs_gather_json(contigs_json)

        # track contigs that have been eliminated at various ranks
        eliminate = set()
        print(f'examining {genome_name} for contamination:')

        summary_d[genome_name] = {}
        summary_d[genome_name]['f_ident'] = f_ident
        summary_d[genome_name]['f_major'] = f_major
        summary_d[genome_name]['comment'] = comment
        summary_d[genome_name]['lineage'] = genome_lineage
        summary_d[genome_name]['filter_at'] = filter_at

        total_bad_n = 0
        total_bad_bp = 0
        for rank in sourmash.lca.taxlist():
            (good_names, bad_names) = calculate_contam(genome_lineage,
                                                       contigs_d,
                                                       rank,
                                                       filter_names=eliminate)
            eliminate.update(bad_names)
            bad_n = len(bad_names)
            bad_bp = sum(bad_names.values())
            total_bad_n += bad_n
            total_bad_bp += bad_bp

            print(f'   {rank}: {len(bad_names)} contigs w/ {kb(bad_bp)}kb')
            summary_d[genome_name][rank] = total_bad_bp
            if rank == match_rank:
                break

        summary_d[genome_name]['total_bad_bp'] = total_bad_bp

        print(f'   (total): {total_bad_n} contigs w/ {kb(total_bad_bp)}kb')

    # output a sorted summary CSV
    fp = open(args.output, 'wt')
    summary_w = csv.writer(fp)
    
    summary_w.writerow(['genome', 'filter_at', 'override_filter_at',
        'total_bad_bp', 'superkingdom_bad_bp', 'phylum_bad_bp',
        'class_bad_bp', 'order_bad_bp', 'family_bad_bp', 'genus_bad_bp',
        'f_ident', 'f_major', 'lineage', 'comment'])

    summary_items = list(summary_d.items())
    summary_items.sort(key=lambda x: -x[1]["total_bad_bp"])

    for genome_name, vals in summary_items:
        vals = summary_d[genome_name]
        summary_w.writerow([genome_name,
                            vals['filter_at'], '',
                            vals["total_bad_bp"],
                            vals['superkingdom'],
                            vals['phylum'],
                            vals['class'],
                            vals['order'],
                            vals['family'],
                            vals['genus'],
                            f'{vals["f_ident"]:.03}',
                            f'{vals["f_major"]:.03}',
                            sourmash.lca.display_lineage(vals["lineage"]),
                            vals["comment"]])

    ####

    print(f"processed {len(genome_names)} genomes.")

    return 0


def cmdline(sys_args):
    "Command line entry point w/argparse action."
    p = argparse.ArgumentParser(sys_args)
    p.add_argument('--input-directory', required=True)
    p.add_argument('--genome-list-file', required=True)
    p.add_argument('-o', '--output', required=True)
    p.add_argument('--lineages-csv', help='lineage spreadsheet', required=True)
    p.add_argument('--provided-lineages', help='provided lineages')
    p.add_argument('--min_f_ident', type=float, default=F_IDENT_THRESHOLD)
    p.add_argument('--min_f_major', type=float, default=F_MAJOR_THRESHOLD)
    args = p.parse_args()

    return main(args)


# execute this, when run with `python -m`.
if __name__ == '__main__':
    returncode = cmdline(sys.argv[1:])
    sys.exit(returncode)
