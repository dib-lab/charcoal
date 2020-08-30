"""
Clean a genome based on hit list & contigs taxonomy.
"""
import sys
import argparse
import gzip
import os.path

import screed

import sourmash

from . import utils
from .version import version
from .utils import (summarize_at_rank, load_contigs_gather_json,
                    is_contig_contaminated, HitList, make_lineage)


def yield_names_in_records(d):
    for k in d:
        r = screed.screedRecord.Record(name=k)
        yield r


def main(args):
    "Main entry point for scripting. Use cmdline for command line entry."

    genome_name = os.path.basename(args.genome)

    hit_list = HitList(args.hit_list)

    try:
        row = hit_list[genome_name]
        filter_at = row['filter_at']
        override_filter_at = row['override_filter_at']
        lineage = row['lineage']
    except KeyError:
        print(f'genome {genome_name} not found in hit list spreadsheet {args.hit_list}; exiting')
        return -1

    filter_rank = filter_at
    if override_filter_at:
        filter_rank = override_filter_at

    if not lineage:
        if filter_rank != 'none':
            print(f'no genome lineage for {genome_name}; cannot clean.')
            print(f"please set filter rank to 'none' rather than {filter_rank}")
            return -1

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
        if not args.do_nothing:
            for record in screed.open(args.genome):
                clean_fp.write(f'>{record.name}\n{record.sequence}\n')
                total_bp += len(record.sequence)
        else:
            total_bp = sum([ x[0] for x in contigs_d.values() ])

        print(f'wrote {total_bp} clean bp to {args.clean}')
        return 0

    print(f'filtering {genome_name} contigs at {filter_rank}')
    print(f'note, genome lineage is {sourmash.lca.display_lineage(lineage)}')

    bp_dirty = 0
    bp_clean = 0

    if not args.do_nothing:
        screed_iter = screed.open(args.genome)
    if args.do_nothing:
        screed_iter = yield_names_in_records(contigs_d)

    for record in screed_iter:
        gather_info = contigs_d[record.name]

        if is_contig_contaminated(lineage, gather_info.gather_tax,
                                  filter_rank, 3): # @CTB configurable?!
            if not args.do_nothing:
                assert len(record.sequence) == gather_info.length
                dirty_fp.write(f'>{record.name}\n{record.sequence}\n')
            bp_dirty += gather_info.length
        else:
            if not args.do_nothing:
                assert len(record.sequence) == gather_info.length
                clean_fp.write(f'>{record.name}\n{record.sequence}\n')
            bp_clean += gather_info.length

    print(f'wrote {bp_clean} clean bp to {args.clean}')
    print(f'wrote {bp_dirty} dirty bp to {args.dirty}')

    return 0


def cmdline(sys_args):
    "Command line entry point w/argparse action."
    p = argparse.ArgumentParser(sys_args)
    p.add_argument('--genome', help='genome file', required=True)
    p.add_argument('--hit-list', help='hit list spreadsheet', required=True)
    p.add_argument('--contigs-json', help='JSON-format contigs classification output by contigs_search', required=True)
    p.add_argument('--clean', help='cleaned contigs', required=True)
    p.add_argument('--dirty', help='dirty contigs', required=True)
    p.add_argument('-n', '--do-nothing', help='do not read or write FASTA')
    args = p.parse_args()

    return main(args)


# execute this, when run with `python -m`.
if __name__ == '__main__':
    returncode = cmdline(sys.argv[1:])
    sys.exit(returncode)
