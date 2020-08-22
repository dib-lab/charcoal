"Tests snakemake execution via click CLI module."
import pytest
import tempfile
import shutil
import os

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


@pytest.mark.dependency()
def test_dory_build_cdbg():
    global _tempdir

    demo_conf = utils.relative_file('demo/demo.conf')
    target = 'output.demo/TOBG_NAT-167.fna.gz.report.txt'
    status = run_snakemake(demo_conf, no_use_conda=True,
                           verbose=True, #outdir=_tempdir,
                           extra_args=[target])
    assert status == 0
    assert os.path.exists(os.path.join(_tempdir, target))


#@pytest.mark.dependency(depends=['test_dory_build_cdbg'])
#def test_dory_build_contigs():
#    global _tempdir

#    dory_conf = utils.relative_file('spacegraphcats/conf/dory-test.yaml')
#    target = 'dory_k21_r1/contigs.fa.gz'
#    status = run_snakemake(dory_conf, verbose=True, outdir=_tempdir,
#                           extra_args=[target])
#    assert status == 0
#    assert os.path.exists(os.path.join(_tempdir, target))
