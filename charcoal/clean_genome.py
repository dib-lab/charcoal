"""
Clean a genome based on hit list & contigs taxonomy.
"""
import sys
import argparse
import gzip
import csv
import os.path

import screed

import sourmash
from sourmash.lca import LineagePair, taxlist

from . import utils
from .version import version
from .utils import (summarize_at_rank, load_contigs_gather_json)


def make_lineage(lineage):
    "Turn a ; or ,-separated set of lineages into a tuple of LineagePair objs."
    lin = lineage.split(';')
    if len(lin) == 1:
        lin = lineage.split(',')
    lin = [ LineagePair(rank, n) for (rank, n) in zip(taxlist(), lin) ]

    return lin


def main(args):
    "Main entry point for scripting. Use cmdline for command line entry."

    genome_name = os.path.basename(args.genome)
    found = False
    with open(args.hit_list, 'rt') as fp:
        r = csv.DictReader(fp)
        for row in r:
            genome = row['genome']
            filter_at = row['filter_at']
            override_filter_at = row['override_filter_at']
            lineage = row['lineage']

            if genome == genome_name:
                found = True
                break

    if not found:
        print(f'genome {genome_name} not found in hit list spreadsheet {args.hit_list}; exiting')
        sys.exit(-1)

    filter_rank = filter_at
    if override_filter_at:
        filter_rank = override_filter_at

    if not lineage:
        if filter_rank != 'none':
            print(f'no genome lineage for {genome_name}; cannot clean.')
            print(f"please set filter rank to 'none' rather than {filter_rank}")
            sys.exit(-1)

    lineage = make_lineage(lineage)

    # load contigs JSON file
    contigs_d = load_contigs_gather_json(args.contigs_json)
    print(f'loaded {len(contigs_d)} contig assignments.')

    xopen = open
    if args.clean.endswith('.gz'):
        xopen = gzip.open
    clean_fp = xopen(args.clean, 'wt')

    xopen = open
    if args.clean.endswith('.gz'):
        xopen = gzip.open
    dirty_fp = xopen(args.dirty, 'wt')

    if filter_rank == 'none':
        dirty_fp.close()

        print(f'filter rank is {filter_rank}; not doing any cleaning.')
        total_bp = 0
        for record in screed.open(args.genome):
            clean_fp.write(f'>{record.name}\n{record.sequence}\n')
            total_bp += len(record.sequence)

        print(f'wrote {total_bp} clean bp to {args.clean}')
        sys.exit(0)

    print(f'filtering {genome_name} contigs at {filter_rank}')
    print(f'note, genome lineage is {sourmash.lca.display_lineage(lineage)}')

    bp_dirty = 0
    bp_clean = 0
    for record in screed.open(args.genome):
        contig_len, contig_taxlist = contigs_d[record.name]
        contig_taxlist_at_rank = summarize_at_rank(contig_taxlist, filter_rank)
        top_hit = None
        if contig_taxlist_at_rank:
            top_hit = contig_taxlist_at_rank[0][0]

        if top_hit and not utils.is_lineage_match(lineage, top_hit, filter_rank):
            dirty_fp.write(f'>{record.name}\n{record.sequence}\n')
            bp_dirty += len(record.sequence)
        else:
            clean_fp.write(f'>{record.name}\n{record.sequence}\n')
            bp_clean += len(record.sequence)

    print(f'wrote {bp_clean} clean bp to {args.clean}')
    print(f'wrote {bp_dirty} dirty bp to {args.dirty}')

def cmdline(sys_args):
    "Command line entry point w/argparse action."
    p = argparse.ArgumentParser(sys_args)
    p.add_argument('--genome', help='genome file', required=True)
    p.add_argument('--hit_list', help='hit list spreadsheet', required=True)
    p.add_argument('--contigs_json', help='JSON-format contigs classification output by contigs_search', required=True)
    p.add_argument('--clean', help='cleaned contigs', required=True)
    p.add_argument('--dirty', help='dirty contigs', required=True)
    args = p.parse_args()

    main(args)


# execute this, when run with `python -m`.
if __name__ == '__main__':
    returncode = cmdline(sys.argv[1:])
    sys.exit(0)
