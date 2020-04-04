#! /usr/bin/env python
"""
Assign taxonomy to shredded fragments in many genomes.

This does the same thing as genome_shred_to_tax, but for many genomes at once.
"""
import sys
import argparse
from pickle import dump
import csv
from collections import defaultdict
import os

import sourmash
import screed
from sourmash.lca import lca_utils
from . import utils                              # charcoal utils
from .genome_shred_to_tax import summarize, classify_signature, shred_to_tax


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--lca-db', required=True)
    p.add_argument('--genomes', nargs='+', required=True)
    p.add_argument('--csv-output-template', default=None)
    p.add_argument('--fragment', default=100000, type=int)
    p.add_argument('--save-tax-hashes-template', default=None)
    args = p.parse_args()

    assert args.csv_output_template

    db, ksize, scaled = lca_utils.load_single_database(args.lca_db)
    mh_factory = sourmash.MinHash(n=0, ksize=ksize, scaled=scaled)
    print('** LCA database:', args.lca_db, ksize, scaled)

    for genome in args.genomes:
        genome_base = os.path.basename(genome)
        output = args.csv_output_template.format(genome=genome_base)

        save_tax_hashes = None
        if args.save_tax_hashes_template:
            save_tax_hashes = args.save_tax_hashes_template.format(genome=genome_base)

        shred_to_tax(genome, output, save_tax_hashes, args.fragment, db,
                     args.lca_db, mh_factory)


    return 0


if __name__ == '__main__':
    sys.exit(main())
