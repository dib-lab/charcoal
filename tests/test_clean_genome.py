import sys
import io
import os.path
from . import pytest_utils as utils

from charcoal import clean_genome


@utils.in_tempdir
def test_1_loomba(location):
    # regression test/check for same results on Loomba
    args = utils.Args()
    args.genome = utils.relative_file("demo/genomes/LoombaR_2017__SID1050_bax__bin.11.fa.gz")
    args.hit_list = utils.relative_file("tests/test-data/loomba-hit-list.csv")
    args.contigs_json = utils.relative_file("tests/test-data/loomba/LoombaR_2017__SID1050_bax__bin.11.fa.gz.contigs-tax.json")
    args.do_nothing = True

    clean_out = os.path.join(location, 'clean.fa.gz')
    dirty_out = os.path.join(location, 'dirty.fa.gz')
    args.clean = clean_out
    args.dirty = dirty_out
    
    status = clean_genome.main(args)

    assert status == 0
    assert os.path.exists(clean_out)
    assert os.path.exists(dirty_out)


@utils.in_tempdir
def test_1_loomba_genus(location):
    # regression test/check for same results on Loomba
    args = utils.Args()
    args.genome = utils.relative_file("demo/genomes/LoombaR_2017__SID1050_bax__bin.11.fa.gz")
    args.hit_list = utils.relative_file("tests/test-data/loomba-hit-list.csv")
    args.contigs_json = utils.relative_file("tests/test-data/loomba/LoombaR_2017__SID1050_bax__bin.11.fa.gz.contigs-tax.json")
    args.do_nothing = True

    clean_out = os.path.join(location, 'clean.fa.gz')
    dirty_out = os.path.join(location, 'dirty.fa.gz')
    args.clean = clean_out
    args.dirty = dirty_out

    sys.stdout, old_out = io.StringIO(), sys.stdout
    sys.stderr, old_err = io.StringIO(), sys.stderr
    try:
        status = clean_genome.main(args)
    finally:
        sys.stdout, new_out = old_out, sys.stdout
        sys.stderr, new_err = old_err, sys.stderr

    assert status == 0
    assert os.path.exists(clean_out)
    assert os.path.exists(dirty_out)

    out = new_out.getvalue()
    assert 'wrote 2727927 clean bp' in out
    assert 'wrote 12351 dirty bp' in out


@utils.in_tempdir
def test_1_loomba_family(location):
    # regression test/check for same results on Loomba
    args = utils.Args()
    args.genome = utils.relative_file("demo/genomes/LoombaR_2017__SID1050_bax__bin.11.fa.gz")

    # swap genus for family in the hit list --
    hit_list = os.path.join(location, 'loomba-hit-list.csv')
    with open(utils.relative_file("tests/test-data/loomba-hit-list.csv"), 'rt') as fp:
        data = fp.read()
        data = data.replace(',genus,', ',family,')

        with open(hit_list, 'wt') as fp2:
            fp2.write(data)

    args.hit_list = hit_list
    args.contigs_json = utils.relative_file("tests/test-data/loomba/LoombaR_2017__SID1050_bax__bin.11.fa.gz.contigs-tax.json")
    args.do_nothing = True

    clean_out = os.path.join(location, 'clean.fa.gz')
    dirty_out = os.path.join(location, 'dirty.fa.gz')
    args.clean = clean_out
    args.dirty = dirty_out

    sys.stdout, old_out = io.StringIO(), sys.stdout
    sys.stderr, old_err = io.StringIO(), sys.stderr
    try:
        status = clean_genome.main(args)
    finally:
        sys.stdout, new_out = old_out, sys.stdout
        sys.stderr, new_err = old_err, sys.stderr

    assert status == 0
    assert os.path.exists(clean_out)
    assert os.path.exists(dirty_out)

    out = new_out.getvalue()
    print(out)
    assert 'wrote 2730931 clean bp' in out
    assert 'wrote 9347 dirty bp' in out


@utils.in_tempdir
def test_1_loomba_order(location):
    # regression test/check for same results on Loomba
    args = utils.Args()
    args.genome = utils.relative_file("demo/genomes/LoombaR_2017__SID1050_bax__bin.11.fa.gz")

    # swap genus for family in the hit list --
    hit_list = os.path.join(location, 'loomba-hit-list.csv')
    with open(utils.relative_file("tests/test-data/loomba-hit-list.csv"), 'rt') as fp:
        data = fp.read()
        data = data.replace(',genus,', ',order,')

        with open(hit_list, 'wt') as fp2:
            fp2.write(data)

    args.hit_list = hit_list
    args.contigs_json = utils.relative_file("tests/test-data/loomba/LoombaR_2017__SID1050_bax__bin.11.fa.gz.contigs-tax.json")
    args.do_nothing = True

    clean_out = os.path.join(location, 'clean.fa.gz')
    dirty_out = os.path.join(location, 'dirty.fa.gz')
    args.clean = clean_out
    args.dirty = dirty_out

    sys.stdout, old_out = io.StringIO(), sys.stdout
    sys.stderr, old_err = io.StringIO(), sys.stderr
    try:
        status = clean_genome.main(args)
    finally:
        sys.stdout, new_out = old_out, sys.stdout
        sys.stderr, new_err = old_err, sys.stderr

    assert status == 0
    assert os.path.exists(clean_out)
    assert os.path.exists(dirty_out)

    out = new_out.getvalue()
    print(out)
    assert 'wrote 2732992 clean bp' in out
    assert 'wrote 7286 dirty bp' in out


@utils.in_tempdir
def test_1_loomba_class(location):
    # regression test/check for same results on Loomba
    args = utils.Args()
    args.genome = utils.relative_file("demo/genomes/LoombaR_2017__SID1050_bax__bin.11.fa.gz")

    # swap genus for family in the hit list --
    hit_list = os.path.join(location, 'loomba-hit-list.csv')
    with open(utils.relative_file("tests/test-data/loomba-hit-list.csv"), 'rt') as fp:
        data = fp.read()
        data = data.replace(',genus,', ',class,')

        with open(hit_list, 'wt') as fp2:
            fp2.write(data)

    args.hit_list = hit_list
    args.contigs_json = utils.relative_file("tests/test-data/loomba/LoombaR_2017__SID1050_bax__bin.11.fa.gz.contigs-tax.json")
    args.do_nothing = True

    clean_out = os.path.join(location, 'clean.fa.gz')
    dirty_out = os.path.join(location, 'dirty.fa.gz')
    args.clean = clean_out
    args.dirty = dirty_out

    sys.stdout, old_out = io.StringIO(), sys.stdout
    sys.stderr, old_err = io.StringIO(), sys.stderr
    try:
        status = clean_genome.main(args)
    finally:
        sys.stdout, new_out = old_out, sys.stdout
        sys.stderr, new_err = old_err, sys.stderr

    print(new_out.getvalue())

    assert status == 0
    assert os.path.exists(clean_out)
    assert os.path.exists(dirty_out)

    out = new_out.getvalue()
    assert 'wrote 2740278 clean bp' in out
    assert 'wrote 0 dirty bp' in out
