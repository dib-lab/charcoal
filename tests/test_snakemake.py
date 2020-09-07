"Tests snakemake execution via click CLI module."
import pytest
import tempfile
import shutil
import os
from pytest_dependency import depends
import io
import sys

from charcoal.__main__ import run_snakemake
from . import pytest_utils as utils

# NOTE re dependencies (@pytest.mark.dependency):
# - These basically duplicate the snakemake dependencies.
# - they're there for convenience, because...
# - ...if these are wrong, the tests will still succeed, they just may
#   do some extra work in some tests & take longer.


def setup_module(m):
    global _tempdir
    _tempdir = tempfile.mkdtemp(prefix='charcoal_test')


def teardown_module(m):
    global _tempdir
    try:
        shutil.rmtree(_tempdir, ignore_errors=True)
    except OSError:
        pass

def _run_snakemake_test(conf, target, extra_args=[]):
    conf = utils.relative_file(conf)
    target = os.path.join(_tempdir, target)

    sys.stdout, old_out = io.StringIO(), sys.stdout
    sys.stderr, old_err = io.StringIO(), sys.stderr
    try:
        status = run_snakemake(conf, no_use_conda=True, verbose=True,
                               outdir=_tempdir,
                               extra_args=[target] + extra_args)
    finally:
        sys.stdout, new_out = old_out, sys.stdout
        sys.stderr, new_err = old_err, sys.stderr

    return status, new_out.getvalue(), new_err.getvalue()


demo_genomes = ['GCF_000005845-subset.fa.gz',
                'LoombaR_2017__SID1050_bax__bin.11.fa.gz',
                'TOBG_NAT-167.fna.gz',
                'TARA_ANE_MAG_00014.fa.gz',
                'TARA_PON_MAG_00084.fa.gz',
                'GCA_001593925.1_ASM159392v1_genomic.fna.gz']


### generate .sig files for each demo genome
@pytest.mark.dependency()
@pytest.mark.parametrize("genome_file", demo_genomes)
def test_make_sig(genome_file):
    target = f'{genome_file}.sig'
    status, out, err = _run_snakemake_test('demo/demo.conf', target)

    assert status == 0
    assert os.path.exists(os.path.join(_tempdir, target))


### generate search matches for each demo genome
@pytest.mark.dependency()
@pytest.mark.parametrize("genome_file", demo_genomes)
def test_make_gather_matches(request, genome_file):
    depends(request, [f"test_make_sig[{g}]" for g in demo_genomes])
    target = f'{genome_file}.matches.sig'
    status, out, err = _run_snakemake_test('demo/demo.conf', target)

    assert status == 0
    assert os.path.exists(os.path.join(_tempdir, target))


### generate JSON contig files for each demo genome
@pytest.mark.dependency()
@pytest.mark.parametrize("genome_file", demo_genomes)
def test_make_contigs_json(request, genome_file):
    # make sure the .sig prereq was generated by test_make_sig!
    depends(request, [f"test_make_gather_matches[{genome_file}]"])

    target = f'{genome_file}.contigs-tax.json'
    status, out, err = _run_snakemake_test('demo/demo.conf', target)
    assert status == 0
    assert os.path.exists(os.path.join(_tempdir, target))


@pytest.mark.dependency()
def test_make_hit_list_dna(request):
    depends(request, [f"test_make_contigs_json[{g}]" for g in demo_genomes])

    target = 'hit_list_for_filtering.csv'
    status, out, err = _run_snakemake_test('demo/demo.conf', target, ['-j', '4'])

    assert status == 0
    assert os.path.exists(os.path.join(_tempdir, target))

@pytest.mark.dependency(['test_make_hit_list_dna'])
@pytest.mark.parametrize("genome_file", demo_genomes)
def test_make_clean_dna(genome_file):
    target = f'{genome_file}.clean.fa.gz'
    status, out, err = _run_snakemake_test('demo/demo.conf', target, ['-j', '4'])

    print(out)
    print(err)

    assert status == 0
    assert os.path.exists(os.path.join(_tempdir, target))


@pytest.mark.dependency()
def test_make_output_prot(request):
#    depends(request, [f"test_make_report[{g}]" for g in demo_genomes])

    target = 'hit_list_for_filtering.csv'
    status, out, err = _run_snakemake_test('demo/demo.prot.conf', target, ['-j', '4'])

    assert status == 0
    assert os.path.exists(os.path.join(_tempdir, target))
