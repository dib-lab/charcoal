import os.path
import shutil
from . import pytest_utils as utils

from charcoal import compare_taxonomy


@utils.in_tempdir
def test_basic(location):
    genome_list_file = os.path.join(location, 'genome-list.txt')
    with open(genome_list_file, 'wt') as fp:
        fp.write("loomba\n")

    output = os.path.join(location, 'hitlist.csv')

    args = utils.Args()
    args.input_directory = location
    args.genome_list_file = genome_list_file
    args.output = output
    args.lineages_csv = utils.relative_file("tests/test-data/test-match-lineages.csv")
    args.provided_lineages = None
    args.min_f_ident = compare_taxonomy.F_IDENT_THRESHOLD
    args.min_f_major = compare_taxonomy.F_MAJOR_THRESHOLD
    
    shutil.copyfile(utils.relative_file("tests/test-data/loomba/LoombaR_2017__SID1050_bax__bin.11.fa.gz.sig"), os.path.join(location, "loomba.sig"))
    shutil.copyfile(utils.relative_file("tests/test-data/loomba/LoombaR_2017__SID1050_bax__bin.11.fa.gz.matches.sig"), os.path.join(location, "loomba.matches.sig"))
    shutil.copyfile(utils.relative_file("tests/test-data/loomba/LoombaR_2017__SID1050_bax__bin.11.fa.gz.contigs-tax.json"), os.path.join(location, "loomba.contigs-tax.json"))

    status = compare_taxonomy.main(args)

    assert status == 0
    assert os.path.exists(output)

    with open(output, 'rt') as fp:
        output_csv = fp.read()

    assert 'loomba,genus,,12351,0,0,0,7286,9347,12351,0.959,0.764,d__Bacteria;p__Firmicutes_A;c__Clostridia;o__Oscillospirales;f__Acutalibacteraceae;g__Anaeromassilibacillus,' in output_csv
