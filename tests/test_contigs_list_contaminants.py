import os.path
from . import pytest_utils as utils
import json

from charcoal import contigs_list_contaminants


@utils.in_tempdir
def test_1_loomba(location):
    # regression test/check for same results on Loomba
    args = utils.Args()
    args.genome = utils.relative_file("demo/genomes/LoombaR_2017__SID1050_bax__bin.11.fa.gz")
    args.genome_sig = utils.relative_file("tests/test-data/loomba/LoombaR_2017__SID1050_bax__bin.11.fa.gz.sig")
    args.matches_csv = utils.relative_file("tests/test-data/loomba/LoombaR_2017__SID1050_bax__bin.11.fa.gz.matches.csv")
    args.databases = [utils.relative_file('tests/test-data/loomba/LoombaR_2017__SID1050_bax__bin.11.fa.gz.matches.zip')]
    args.lineages_csv = utils.relative_file("tests/test-data/test-match-lineages.csv")
    args.hitlist = utils.relative_file("tests/test-data/loomba-hit-list.csv")
    args.json_out = os.path.join(location, 'tax.json')
    args.match_rank = 'genus'

    status = contigs_list_contaminants.main(args)

    assert status == 0
    assert os.path.exists(args.json_out)

    with open(args.json_out, 'rt') as fp:
        results = json.load(fp)
        assert results != {}


@utils.in_tempdir
def test_1_loomba_abund(location):
    # regression test/check for same results on Loomba - with abund sigs
    args = utils.Args()
    args.genome = utils.relative_file("demo/genomes/LoombaR_2017__SID1050_bax__bin.11.fa.gz")
    args.genome_sig = utils.relative_file("tests/test-data/loomba/LoombaR_2017__SID1050_bax__bin.11.fa.gz.sig")
    args.matches_csv = utils.relative_file("tests/test-data/loomba/LoombaR_2017__SID1050_bax__bin.11.fa.gz.matches.csv")
    args.databases = [utils.relative_file('tests/test-data/loomba/LoombaR_2017__SID1050_bax__bin.11.fa.gz.matches.abund.zip')]
    args.lineages_csv = utils.relative_file("tests/test-data/test-match-lineages.csv")
    args.hitlist = utils.relative_file("tests/test-data/loomba-hit-list.csv")
    args.json_out = os.path.join(location, 'tax.json')
    args.match_rank = 'genus'

    status = contigs_list_contaminants.main(args)

    assert status == 0
    assert os.path.exists(args.json_out)

    with open(args.json_out, 'rt') as fp:
        results = json.load(fp)
        assert results != {}
