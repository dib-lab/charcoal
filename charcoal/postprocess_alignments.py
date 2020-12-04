#! /usr/bin/env python
"""
Create a "hit list" of how much will be removed at what ranks.
"""
import sys
import argparse
import csv
import os.path
import yaml

import sourmash
from . import utils

def main(args):
    "Main entry point for scripting. Use cmdline for command line entry."

    inp_dir = args.input_directory
    hitlist = utils.CSV_DictHelper(args.hit_list, 'genome')

    with open(args.matches_yaml, 'rt') as fp:
        matches_info = yaml.safe_load(fp)

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
