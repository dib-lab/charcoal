#! /usr/bin/env python
"""
Create a "hit list" of how much will be removed at what ranks.
"""
import sys
import argparse
import csv
import os.path
from collections import Counter
import json

import sourmash
from sourmash.lca.command_index import load_taxonomy_assignments
from sourmash.lca import LCA_Database, LineagePair

from . import utils
from .lineage_db import LineageDB
from .utils import (gather_at_rank, get_ident, summarize_at_rank,
                    pretty_print_lineage, load_contigs_gather_json,
                    is_contig_contaminated, is_contig_clean)


GATHER_MIN_MATCHES=3
F_IDENT_THRESHOLD=0.1
F_MAJOR_THRESHOLD=0.2


def kb(bp):
    return int(bp/1000)


def calculate_contam(genome_lin, contigs_d, rank, filter_names=None):
    "Calculate not-bad bp at each rank. Be conservative."
    good_names = dict()
    bad_names = dict()

    for contig_name, gather_info in contigs_d.items():
        contig_taxlist = gather_info.gather_tax
        if filter_names and contig_name in filter_names:
            continue

        if is_contig_contaminated(genome_lin, contig_taxlist, rank, GATHER_MIN_MATCHES):
            bad_names[contig_name] = gather_info
        else:
            good_names[contig_name] = gather_info

    return (good_names, bad_names)


def calculate_clean(genome_lin, contigs_d, rank):
    "Calculate definitely-clean bp, as opposed to not-bad bp."
    good_names = dict()
    bad_names = dict()

    for contig_name, gather_info in contigs_d.items():
        contig_taxlist = gather_info.gather_tax

        if not is_contig_contaminated(genome_lin, contig_taxlist, rank, GATHER_MIN_MATCHES):
            good_names[contig_name] = gather_info
        else:
            bad_names[contig_name] = gather_info

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
    
    f_ident = 0
    f_major = 0
    if sum_ident:
        f_ident = sum_ident / len(entire_mh)
        f_major = first_count / sum_ident

    return first_lin, f_ident, f_major


def choose_genome_lineage(guessed_genome_lineage, provided_lineage, match_rank,
                          f_ident, f_major, min_f_ident, min_f_major):

    comment = ""
    genome_lineage = None
    needs_lineage = False

    if provided_lineage:
        if utils.is_lineage_match(provided_lineage, guessed_genome_lineage, match_rank):
            print(f'(provided lineage agrees with k-mer classification at {match_rank} level)')
        elif guessed_genome_lineage:
            print(f'(provided lineage disagrees with k-mer classification at or above {match_rank} level)')
        else:
            pass

        genome_lineage = utils.pop_to_rank(provided_lineage, match_rank)
        print(f'\nUsing provided lineage as genome lineage.')
    else:
        if f_ident < min_f_ident:
            print(f'** ERROR: fraction of total identified hashes (f_ident) < {min_f_ident*100:.0f}%.')
            comment = f"too few identifiable hashes; f_ident < {min_f_ident*100:.0f}%. provide a lineage for this genome."
            needs_lineage = True
        elif f_major < min_f_major:
            print(f'** ERROR: fraction of identified hashes in major lineage (f_major) < {min_f_major*100:.0f}%.')
            comment = f"too few hashes in major lineage; f_major < {min_f_major*100:.0f}%. provide a lineage for this genome."
            needs_lineage = True
        else:
            genome_lineage = utils.pop_to_rank(guessed_genome_lineage, match_rank)
            print(f'Using majority gather lineage as genome lineage.')

    return genome_lineage, comment, needs_lineage


def get_genome_taxonomy(matches_filename, database_list,
                        genome_sig_filename, provided_lineage,
                        tax_assign, match_rank, min_f_ident, min_f_major):
    # load the matches from prefetch as a picklist
    picklist = sourmash.picklist.SignaturePicklist('prefetch')
    try:
        picklist.load(matches_filename, picklist.column_name)
    except ValueError:
        with open(matches_filename, 'rt') as fp:
            contents = fp.read()
            if not len(contents): # empty is ok.
                picklist = None
            else:
                raise

    # load all of the matches in the database, as found by prefetch;
    # select on them; and then aggregate into MultiIndex.
    # CTB note: currently, this loads all the signatures into memory.
    # Alternatively we could do something with LazyLoadedIndex maybe?

    siglist = []
    for filename in database_list:
        db = sourmash.load_file_as_index(filename)
        db = db.select(picklist=picklist)
        siglist += list(db.signatures())

    if not siglist:
        comment = 'no matches for this genome.'
        print(comment)
        return None, comment, False, 0.0, 0.0

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
                print(f'found exact match: {ss.name}. removing.')
            else:
                print(f'found exact match: {ss.name}. but no provided lineage!')
                comment = f'found exact match: {ss.name}. but no provided lineage! cannot analyze.'
                return None, comment, True, 1.0, 1.0

    # ...but leave exact matches in if they're the only matches, I guess!
    if new_siglist:
        siglist = new_siglist

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
        print(comment)
        return guessed_genome_lineage, comment, False, f_major, f_ident

    # did we get a passed-in lineage assignment?
    provided_lin = ""
    if provided_lineage:
        provided_lin = [ LineagePair(rank, name) for (rank, name) in zip(sourmash.lca.taxlist(), provided_lineage) if name.strip() ]
        print(f'Provided lineage from command line:\n   {sourmash.lca.display_lineage(provided_lin)}')

    # choose between the lineages
    genome_lineage, comment, needs_lineage = \
                              choose_genome_lineage(guessed_genome_lineage,
                                                    provided_lin,
                                                    match_rank,
                                                    f_ident, f_major,
                                                    min_f_ident, min_f_major)

    return genome_lineage, comment, needs_lineage, f_major, f_ident

###

def main(args):
    "Main entry point for scripting. Use cmdline for command line entry."
    match_rank = args.match_rank

    # load taxonomy assignments for all the things
    tax_assign, _ = load_taxonomy_assignments(args.lineages_csv,
                                              start_column=2)
    print(f'loaded {len(tax_assign)} tax assignments.')

    # place to load in the genome from
    dirname = args.input_directory

    genome_name = args.genome
    print(f"working on {genome_name}")

    # load the provided lineages file
    provided_lineages = {}
    if args.provided_lineages:
        with open(args.provided_lineages, 'rt') as fp:
            r = csv.reader(fp)
            for row in r:
                gname = row[0]
                lin = row[1:]
                provided_lineages[gname] = lin

    print(f"loaded {len(provided_lineages)} provided lineages")

    ### XXX CTB

    detected_contam = {}

    summary_d = {}
    matches_filename = os.path.join(dirname, genome_name + '.matches.csv')
    genome_sig = os.path.join(dirname, genome_name + '.sig')
    lineage = provided_lineages.get(genome_name, '')
    contigs_json = os.path.join(dirname, genome_name + '.contigs-tax.json')

    x = get_genome_taxonomy(matches_filename,
                            args.databases,
                            genome_sig,
                            lineage,
                            tax_assign, match_rank,
                            args.min_f_ident,
                            args.min_f_major)
    genome_lineage, comment, needs_lineage, f_major, f_ident = x

    # did we get a lineage for this genome? if so, propose filtering at
    # default rank 'match_rank', otherwise ...do not filter.
    if genome_lineage and not comment:
        filter_at = match_rank
    else:
        genome_lineage = []
        filter_at = 'none'

    # load contigs tax
    contigs_d = load_contigs_gather_json(contigs_json)

    print(f'examining {genome_name} for contamination')

    vals = {}
    vals['genome'] = genome_name
    vals['f_ident'] = f_ident
    vals['f_major'] = f_major
    vals['comment'] = comment
    vals['needs_lineage_flag'] = 1 if needs_lineage else 0
    vals['lineage'] = sourmash.lca.display_lineage(genome_lineage)
    vals['filter_at'] = filter_at

    # calculate summary stats for contigs
    nohash_bp = 0
    nohash_count = 0
    noident_bp = 0
    noident_count = 0
    contigs_n = 0
    contigs_bp = 0
    for contig_name, gather_info in contigs_d.items():
        contigs_n += 1
        contigs_bp += gather_info.length

        if not gather_info.num_hashes:
            nohash_bp += gather_info.length
            nohash_count += 1
        elif not gather_info.gather_tax or not genome_lineage:
            noident_bp += gather_info.length
            noident_count += 1

    vals['ignored_contigs_n'] = nohash_count
    vals['ignored_contigs_bp'] = nohash_bp
    vals['noident_contigs_n'] = noident_count
    vals['noident_contigs_bp'] = noident_bp
    vals['total_contigs_n'] = contigs_n
    vals['total_contigs_bp'] = contigs_bp

    # track contigs that have been eliminated at various ranks
    for rank in sourmash.lca.taxlist():
        (good_names, bad_names) = calculate_contam(genome_lineage,
                                                   contigs_d,
                                                   rank)

        bad_n = len(bad_names)
        bad_bp = sum([ x.length for x in bad_names.values() ])

        print(f'   {rank}: {len(bad_names)} contigs w/ {kb(bad_bp)}kb')
        vals[f'bad_{rank}_bp'] = bad_bp
        vals[f'bad_{rank}_n'] = bad_n

        (good_names, bad_names) = calculate_clean(genome_lineage,
                                                  contigs_d,
                                                  rank)
        good_n = len(good_names)
        good_bp = sum([ x.length for x in good_names.values() ])
        vals[f'good_{rank}_bp'] = good_bp
        vals[f'good_{rank}_n'] = good_n

        assert bad_bp + good_bp == contigs_bp

        # track contamination between source (genome) / target (contig)
        x = []
        for contig_name in bad_names:
            contig_taxlist = contigs_d[contig_name].gather_tax
            for hit, count in contig_taxlist:
                if utils.is_lineage_match(genome_lineage, hit, rank):
                    continue

                # contam!
                source_lin = utils.pop_to_rank(genome_lineage, rank)
                target_lin = utils.pop_to_rank(hit, rank)

                x.append((source_lin, target_lin, count))

        detected_contam[genome_name] = x


        if rank == 'genus':
            break

    vals['total_bad_bp'] = vals[f'bad_{match_rank}_bp']
    if vals['total_bad_bp'] == 0:
        vals['filter_at'] = 'none'

    print(f"   (total): {vals['bad_genus_n']} contigs w/ {kb(vals['bad_genus_bp'])}kb")

    summary_d[genome_name] = vals

    # output a hit list CSV for this genome
    fp = open(args.hit_list, 'wt')
    hitlist_w = csv.writer(fp)
    
    hitlist_w.writerow(['genome', 'filter_at', 'override_filter_at',
        'total_bad_bp', 'superkingdom_bad_bp', 'phylum_bad_bp',
        'class_bad_bp', 'order_bad_bp', 'family_bad_bp', 'genus_bad_bp',
        'f_ident', 'f_major', 'lineage', 'comment'])

    summary_items = list(summary_d.items())
    assert len(summary_items) == 1
    #summary_items.sort(key=lambda x: -x[1]["total_bad_bp"])

    for genome_name, vals in summary_items:
        hitlist_w.writerow([genome_name,
                            vals['filter_at'], '',
                            vals["total_bad_bp"],
                            vals['bad_superkingdom_bp'],
                            vals['bad_phylum_bp'],
                            vals['bad_class_bp'],
                            vals['bad_order_bp'],
                            vals['bad_family_bp'],
                            vals['bad_genus_bp'],
                            f'{vals["f_ident"]:.03f}',
                            f'{vals["f_major"]:.03f}',
                            vals["lineage"],
                            vals["comment"]])

    fp.close()

    # output a single-line summary CSV with a lot more information!
    fp = open(args.genome_summary, 'wt')

    # build column list; put genome first
    vals = summary_items[0][1]
    all_columns = set(vals.keys())
    all_columns.remove('genome')
    all_columns = list(sorted(all_columns))
    all_columns = ['genome'] + all_columns

    w = csv.DictWriter(fp, fieldnames=all_columns)
    w.writeheader()

    for genome, vals in summary_items:
        w.writerow(vals)

    fp.close()

    ####

    print(f"processed {genome_name}.")

    print(f"saving contamination summary to {args.contam_summary_json}")
    with open(args.contam_summary_json, 'wt') as fp:
        utils.save_contamination_summary(detected_contam, fp)

    return 0


def cmdline(sys_args):
    "Command line entry point w/argparse action."
    p = argparse.ArgumentParser(sys_args)
    p.add_argument('--input-directory', required=True)
    p.add_argument('--hit-list', required=True)
    p.add_argument('--genome-summary', required=True)
    p.add_argument('--contam-summary-json', required=True)   # @CTB remove?
    p.add_argument('--lineages-csv', help='lineage spreadsheet', required=True)
    p.add_argument('--provided-lineages', help='provided lineages')
    p.add_argument('--min_f_ident', type=float, default=F_IDENT_THRESHOLD)
    p.add_argument('--min_f_major', type=float, default=F_MAJOR_THRESHOLD)
    p.add_argument('--match-rank', required=True)
    p.add_argument('--databases', help='sourmash databases', required=True,
                   nargs='+')
    p.add_argument('genome')
    args = p.parse_args()

    return main(args)


# execute this, when run with `python -m`.
if __name__ == '__main__':
    returncode = cmdline(sys.argv[1:])
    sys.exit(returncode)
