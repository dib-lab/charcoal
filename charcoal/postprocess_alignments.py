#! /usr/bin/env python
"""
Create a "hit list" of how much will be removed at what ranks.
"""
import sys
import argparse
import csv
import os.path
import yaml
import glob
import itertools

from . import alignplot
from .alignplot import AlignmentContainer

import sourmash
from . import utils

# minimum alignment size, in kb
# in practice, this is also bounded by the aligner used - mashmap generally
# won't go below 5kb, for example.
MIN_ALIGN_SIZE=0.5


def main(args):
    "Main entry point for scripting. Use cmdline for command line entry."

    inp_dir = args.input_directory
    hitlist = utils.CSV_DictHelper(args.hit_list, 'genome')

    genomebase = os.path.basename(args.genome)

    with open(args.matches_yaml, 'rt') as fp:
        matches_info = yaml.safe_load(fp)

    ##
    genome_lin = utils.make_lineage(matches_info['query_info']['genome_lineage'])
    match_rank = matches_info['query_info']['match_rank']
    scaled = matches_info['query_info']['scaled']

    clean_accs = []
    clean_accs_d = {}
    dirty_accs = []
    dirty_accs_d = {}
    for match_acc, acc_info in matches_info['matches'].items():
        match_counts = acc_info['counts']
        match_type = acc_info['match_type']
        match_lineage = acc_info['lineage']

        assert match_acc not in clean_accs_d
        assert match_acc not in dirty_accs_d

        if match_type == 'clean':
            clean_accs.append((match_acc, match_lineage, match_counts))
            clean_accs_d[match_acc] = (match_lineage, match_counts)
        elif match_type == 'dirty':
            dirty_accs.append((match_acc, match_lineage, match_counts))
            dirty_accs_d[match_acc] = (match_lineage, match_counts)

    clean_accs.sort(key=lambda x: -x[2])
    dirty_accs.sort(key=lambda x: -x[2])

    output = []

    output.append(f'loaded {len(clean_accs)} clean accs and {len(dirty_accs)} dirty accs')
    output.append('')
    output.append(f'query genome lineage: `{utils.display_lineage(genome_lin)}`\n')

    output.append(f'genomes that match the lineage at {match_rank}:')
    for (match_acc, match_lineage, match_counts) in clean_accs:
        output.append(f'* `{match_acc}` with est {match_counts*scaled} kb;\n`{match_lineage}`')

    output.append('')
    output.append('genomes that do NOT match the lineage:')
    for (match_acc, match_lineage, match_counts) in dirty_accs:
        output.append(f'* `{match_acc}` with est {match_counts*scaled} kb;\n`{match_lineage}`')

    print("\n".join(output))

    def load_target_pairs(match_list):
        pairs = []
        for acc, _, _ in match_list:
            filename = glob.glob(f'genbank_genomes/{acc}*.fna.gz')
            #assert len(filename) == 1, filename # @CTB
            filename = filename[0]
            pairs.append((acc, filename))

        return pairs

    contaminant_pairs = load_target_pairs(dirty_accs)
    clean_pairs = load_target_pairs(clean_accs)

    contigs_by_acc = {}
    contigs_to_acc = {}
    all_sizes = {}
    all_sizes.update(alignplot.load_contig_sizes(args.genome))
    for acc, _, _ in itertools.chain(clean_accs, dirty_accs):
        filename = glob.glob(f'genbank_genomes/{acc}*.fna.gz')
        filename = filename[0]
        sizes = alignplot.load_contig_sizes(filename)
        all_sizes.update(sizes)
        contigs_by_acc[acc] = sizes
        for contig_name in sizes:
            assert contig_name not in contigs_to_acc
            contigs_to_acc[contig_name] = acc

    dirty_alignment = AlignmentContainer(genomebase, args.genome, contaminant_pairs, f'{inp_dir}/hitlist-accessions.info.csv')

    results = {}
    for t_acc, _ in contaminant_pairs:
        mashmap_file = f'{inp_dir}/{genomebase}.x.{t_acc}.mashmap.align'
        results[t_acc] = dirty_alignment._read_mashmap(mashmap_file)
    dirty_alignment.results = results

    print(f'filtering dirty alignments to query size >= 500 and identity >= {args.min_align_pident}%')
    dirty_alignment.filter(query_size=MIN_ALIGN_SIZE, pident=args.min_align_pident)
    print(f'filtering dirty alignments to min query coverage {args.min_query_coverage}')
    dirty_alignment = dirty_alignment.filter_by_query_coverage(args.min_query_coverage)

    sum_dirty_kb = sum(dirty_alignment.calc_shared().values())
    print(f'query genome lineage: {utils.display_lineage(genome_lin)}')
    print(f'**dirty bases: {sum_dirty_kb:.1f}kb of alignments to query genome, across all targets.**')

    all_regions = []
    for t_acc, vv in dirty_alignment.results.items():
        all_regions.extend(vv)
    regions_by_query = alignplot.group_regions_by(all_regions, 'query')
    query_shared = dirty_alignment.calc_shared()

    sum_to_remove = 0
    for k, covered_bases in query_shared.items():
        sum_to_remove += all_sizes[k]
        print(f'\nremoving {all_sizes[k]:.0f}kb with {covered_bases:.0f}kb dirty, contig name {k}.')
        for region in regions_by_query[k]:
            source_acc = contigs_to_acc[region.target]
            source_lin = utils.make_lineage(dirty_accs_d[source_acc][0])
            query_aligned = alignplot.region_size(region, 'query')
            print(f'   {query_aligned:.0f}kb aligns to {source_acc}:{region.target} at {region.pident:.1f}%')
            print(f'   ({utils.display_lineage(source_lin)})')
            disagree_rank = utils.find_disagree_rank(genome_lin, source_lin)
            query_at_rank = utils.pop_to_rank(genome_lin, disagree_rank)[-1].name
            source_at_rank = utils.pop_to_rank(source_lin, disagree_rank)[-1].name
            print(f"   ** disagreement at rank '{disagree_rank}'; genome {query_at_rank}, source {source_at_rank}")

    print(f'\nremoving {sum_to_remove:.0f}kb total in contigs >= {args.min_query_coverage*100:.0f}% dirty, based on alignments at {args.min_align_pident:.0f}% identity over {MIN_ALIGN_SIZE:.2f}kb or more')

    ###

    with open(args.yaml_out, 'wt') as fp:
        pass
    with open(args.summary_csv, 'wt') as fp:
        pass

    return 0


def cmdline(sys_args):
    "Command line entry point w/argparse action."
    p = argparse.ArgumentParser(sys_args)
    p.add_argument('--input-directory', required=True)
    p.add_argument('--hit-list', required=True)
    p.add_argument('--matches-yaml', required=True)
    p.add_argument('--yaml-out', required=True)
    p.add_argument('--summary-csv', required=True)
    p.add_argument('--min-query-coverage', type=float, required=True)
    p.add_argument('--min-align-pident', type=float, required=True)
    p.add_argument('genome')
    args = p.parse_args()

    return main(args)


# execute this, when run with `python -m`.
if __name__ == '__main__':
    returncode = cmdline(sys.argv[1:])
    sys.exit(returncode)
