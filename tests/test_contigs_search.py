import os.path
from . import pytest_utils as utils
import json

from charcoal import contigs_search

@utils.in_tempdir
def test_1(location):
    # test an empty set of matches (once self is removed)
    args = utils.Args()
    args.genome = utils.relative_file("tests/test-data/genomes/2.fa.gz")
    args.genome_sig = utils.relative_file("tests/test-data/genomes/2.fa.gz.sig")
    args.matches_sig = utils.relative_file("tests/test-data/2.fa.gz.gather-matches.sig.gz")
    args.lineages_csv = utils.relative_file("tests/test-data/test-match-lineages.csv")
    args.json_out = os.path.join(location, 'tax.json')

    status = contigs_search.main(args)

    assert status == 0
    assert os.path.exists(args.json_out)

    with open(args.json_out, 'rt') as fp:
        results = json.load(fp)
        assert results == {}
