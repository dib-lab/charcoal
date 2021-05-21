import os.path
from . import pytest_utils as utils
import json

from charcoal import contigs_search_taxonomy


@utils.in_tempdir
def test_1(location):
    # test an empty set of matches (once self is removed)
    args = utils.Args()
    args.genome = utils.relative_file("tests/test-data/genomes/2.fa.gz")
    args.genome_sig = utils.relative_file("tests/test-data/genomes/2.fa.gz.sig")
    args.matches_sig = utils.relative_file("tests/test-data/2.fa.gz.gather-matches.zip")
    args.lineages_csv = utils.relative_file("tests/test-data/test-match-lineages.csv")
    args.json_out = os.path.join(location, 'tax.json')
    args.match_rank = 'genus'

    status = contigs_search_taxonomy.main(args)

    assert status == 0
    assert os.path.exists(args.json_out)

    with open(args.json_out, 'rt') as fp:
        results = json.load(fp)
        assert results == {}


@utils.in_tempdir
def test_2_loomba(location):
    # regression test/check for same results on Loomba
    args = utils.Args()
    args.genome = utils.relative_file("demo/genomes/LoombaR_2017__SID1050_bax__bin.11.fa.gz")
    args.genome_sig = utils.relative_file("tests/test-data/loomba/LoombaR_2017__SID1050_bax__bin.11.fa.gz.sig")
    args.matches_sig = utils.relative_file("tests/test-data/loomba/LoombaR_2017__SID1050_bax__bin.11.fa.gz.matches.zip")
    args.lineages_csv = utils.relative_file("tests/test-data/test-match-lineages.csv")
    args.json_out = os.path.join(location, 'tax.json')
    args.match_rank = 'genus'

    status = contigs_search_taxonomy.main(args)

    assert status == 0
    assert os.path.exists(args.json_out)

    import shutil
    shutil.copyfile(args.json_out, '/tmp/xyz.json')
    with open(args.json_out, 'rt') as fp:
        this_results = json.load(fp)

    saved_results_file = utils.relative_file("tests/test-data/loomba/LoombaR_2017__SID1050_bax__bin.11.fa.gz.contigs-tax.json")
    with open(saved_results_file, 'rt') as fp:
        saved_results = json.load(fp)

    assert this_results == saved_results
