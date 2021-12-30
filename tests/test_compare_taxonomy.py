import os.path
import shutil
from . import pytest_utils as utils

from charcoal import compare_taxonomy
from charcoal.utils import load_contamination_summary


@utils.in_tempdir
def test_basic(location):
    genome_list_file = os.path.join(location, 'genome-list.txt')
    with open(genome_list_file, 'wt') as fp:
        fp.write("loomba\n")

    hitlist = os.path.join(location, 'hitlist.csv')
    summary_csv = os.path.join(location, 'summary.csv')
    contam_json = os.path.join(location, 'contam.json')

    args = utils.Args()
    args.input_directory = location
    args.genome_list_file = genome_list_file
    args.hit_list = hitlist
    args.genome_summary = summary_csv
    args.lineages_csv = utils.relative_file("tests/test-data/test-match-lineages.csv")
    args.databases = [utils.relative_file('tests/test-data/loomba/LoombaR_2017__SID1050_bax__bin.11.fa.gz.matches.zip')]
    args.provided_lineages = None
    args.contam_summary_json = contam_json
    args.min_f_ident = compare_taxonomy.F_IDENT_THRESHOLD
    args.min_f_major = compare_taxonomy.F_MAJOR_THRESHOLD
    args.match_rank = 'genus'
    args.genome = 'loomba'
    
    shutil.copyfile(utils.relative_file("tests/test-data/loomba/LoombaR_2017__SID1050_bax__bin.11.fa.gz.sig"), os.path.join(location, "loomba.sig"))
    shutil.copyfile(utils.relative_file("tests/test-data/loomba/LoombaR_2017__SID1050_bax__bin.11.fa.gz.matches.csv"), os.path.join(location, "loomba.matches.csv"))
    shutil.copyfile(utils.relative_file("tests/test-data/loomba/LoombaR_2017__SID1050_bax__bin.11.fa.gz.contigs-tax.json"), os.path.join(location, "loomba.contigs-tax.json"))

    status = compare_taxonomy.main(args)

    assert status == 0
    assert os.path.exists(hitlist)
    assert os.path.exists(contam_json)

    with open(hitlist, 'rt') as fp:
        hitlist_csv = fp.read()

    assert 'loomba,genus,,12351,0,0,0,7286,9347,12351,0.959,0.764,d__Bacteria;p__Firmicutes_A;c__Clostridia;o__Oscillospirales;f__Acutalibacteraceae;g__Anaeromassilibacillus,' in hitlist_csv

    # this next bit is a bit off-topic for this test, but I want to
    # make sure that the saved hit-list matches used in clean_genome
    # tests are kept up to date with the ones actually generated by
    # compare taxonomy, so ... I'm putting this check in the same test :)
    with open(utils.relative_file('tests/test-data/loomba-hit-list.csv'), 'rt') as fp:
        saved_csv = fp.read()

    assert ',genus,,12351,0,0,0,7286,9347,12351,0.959,0.764,d__Bacteria;p__Firmicutes_A;c__Clostridia;o__Oscillospirales;f__Acutalibacteraceae;g__Anaeromassilibacillus,' in saved_csv

    # check contam output against saved
    test_contam_file = utils.relative_file("tests/test-data/loomba/contam.json")
    with open(test_contam_file, 'rt') as fp:
        test_contam = load_contamination_summary(fp)

    with open(args.contam_summary_json, 'rt') as fp:
        actual_contam = load_contamination_summary(fp)

    assert test_contam == actual_contam
