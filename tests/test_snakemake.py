"Tests snakemake execution via click CLI module."
import pytest
import tempfile
import shutil
import os
from pytest_dependency import depends

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

    status = run_snakemake(conf, no_use_conda=True, verbose=True,
                           outdir=_tempdir, extra_args=[target] + extra_args)
    return status


demo_genomes = ['GCF_000005845-subset.fa.gz',
                'LoombaR_2017__SID1050_bax__bin.11.fa.gz',
                'TOBG_NAT-167.fna.gz',
                'TARA_ANE_MAG_00014.fa.gz',
                'TARA_PON_MAG_00084.fa.gz',
                'GCA_001593925.1_ASM159392v1_genomic.fna.gz']


### generate .sig files for each demo genome
### DISABLED for now
@pytest.mark.dependency()
@pytest.mark.parametrize("genome_file", demo_genomes)
def xxx_test_make_sig(genome_file):
    target = f'{genome_file}.sig'
    status = _run_snakemake_test('demo/demo.conf', target)

    assert status == 0
    assert os.path.exists(os.path.join(_tempdir, target))


### generate report.txt files for each demo genome
### DISABLED for now
@pytest.mark.dependency()
@pytest.mark.parametrize("genome_file", demo_genomes)
def xxx_test_make_report(request, genome_file):
    # make sure the .sig prereq was generated by test_make_sig!
    depends(request, [f"test_make_sig[{genome_file}]"])

    target = f'{genome_file}.report.txt'
    status = _run_snakemake_test('demo/demo.conf', target)
    assert status == 0
    assert os.path.exists(os.path.join(_tempdir, target))


@pytest.mark.dependency()
def test_make_output(request):
#    depends(request, [f"test_make_report[{g}]" for g in demo_genomes])

    target = 'combined_summary.csv'
    status = _run_snakemake_test('demo/demo.conf', target, ['-j', '4'])

    assert status == 0
    assert os.path.exists(os.path.join(_tempdir, target))
