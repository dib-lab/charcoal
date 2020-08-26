import os.path
from . import pytest_utils as utils
import json

from charcoal import clean_genome


@utils.in_tempdir
def test_1_loomba(location):
    # regression test/check for same results on Loomba
    args = utils.Args()
    args.genome = utils.relative_file("demo/genomes/LoombaR_2017__SID1050_bax__bin.11.fa.gz")
    args.hit_list = utils.relative_file("tests/test-data/loomba-hit-list.csv")
    args.contigs_json = utils.relative_file("tests/test-data/loomba/LoombaR_2017__SID1050_bax__bin.11.fa.gz.contigs-tax.json")

    clean_out = os.path.join(location, 'clean.fa.gz')
    dirty_out = os.path.join(location, 'dirty.fa.gz')
    args.clean = clean_out
    args.dirty = dirty_out
    
    status = clean_genome.main(args)

    assert status == 0
    assert os.path.exists(clean_out)
    assert os.path.exists(dirty_out)
