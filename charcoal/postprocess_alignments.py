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

from .alignplot import AlignmentContainer

import sourmash
from . import utils

def main(args):
    "Main entry point for scripting. Use cmdline for command line entry."

    inp_dir = args.input_directory
    hitlist = utils.CSV_DictHelper(args.hit_list, 'genome')

    genomebase = os.path.basename(args.genome)

    with open(args.matches_yaml, 'rt') as fp:
        matches_info = yaml.safe_load(fp)

    ##
    genome_lin = matches_info['query_info']['genome_lineage']
    match_rank = matches_info['query_info']['match_rank']
    scaled = matches_info['query_info']['scaled']

    clean_accs = []
    dirty_accs = []
    for match_acc, acc_info in matches_info['matches'].items():
        match_counts = acc_info['counts']
        match_type = acc_info['match_type']
        match_lineage = acc_info['lineage']

        if match_type == 'clean':
            clean_accs.append((match_acc, match_lineage, match_counts))
        elif match_type == 'dirty':
            dirty_accs.append((match_acc, match_lineage, match_counts))

    clean_accs.sort(key=lambda x: -x[2])
    dirty_accs.sort(key=lambda x: -x[2])

    output = []

    output.append(f'loaded {len(clean_accs)} clean accs and {len(dirty_accs)} dirty accs')
    output.append('')
    output.append(f'query genome lineage: `{genome_lin}`\n')

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

    dirty_alignment = AlignmentContainer(genomebase, args.genome, contaminant_pairs, f'{inp_dir}/hitlist-accessions.info.csv')

    results = {}
    for t_acc, _ in contaminant_pairs:
        mashmap_file = f'{inp_dir}/{genomebase}.x.{t_acc}.mashmap.align'
        results[t_acc] = dirty_alignment._read_mashmap(mashmap_file)
    dirty_alignment.results = results

    print('filtering dirty alignments to query size >= 500 and identity >= 95%')
    dirty_alignment.filter(query_size=0.5, pident=95)

    sum_dirty_kb = sum(dirty_alignment.calc_shared().values())
    print(f'**dirty bases: {sum_dirty_kb:.1f}kb of alignments to query genome, across all targets.**')
    

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
    p.add_argument('genome')
    args = p.parse_args()

    return main(args)


# execute this, when run with `python -m`.
if __name__ == '__main__':
    returncode = cmdline(sys.argv[1:])
    sys.exit(returncode)
