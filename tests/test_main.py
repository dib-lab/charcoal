import types
import tempfile
import shutil
import os.path


from charcoal.just_taxonomy import main


def get_test_data(filename):
    filepath = os.path.join(os.path.dirname(__file__), 'test-data',
                            filename)
    return filepath


class TempDirectory(object):
    def __init__(self):
        self.tempdir = tempfile.mkdtemp(prefix='charcoal_test_')

    def __enter__(self):
        return self.tempdir

    def __exit__(self, exc_type, exc_value, traceback):
        try:
            shutil.rmtree(self.tempdir, ignore_errors=True)
        except OSError:
            pass

        if exc_type:
            return False


def build_args(**kwargs):
    args = types.SimpleNamespace(**kwargs)

    # currently always need args.genome, args.match_rank,
    # args.lineages_csv, args.matches_sig, args.summary,
    # args.report, args.contig_report, args.clean, and args.dirty.

    return args


def test_main_1():
    # basic test - run it!
    genome = get_test_data('genomes/2.fa.gz')
    matches_sig = get_test_data('2.fa.gz.gather-matches.sig.gz')
    lineages_csv = get_test_data('test-match-lineages.csv')
    
    with TempDirectory() as tmpdir:
        clean = os.path.join(tmpdir, 'clean.fa.gz')
        dirty = os.path.join(tmpdir, 'dirty.fa.gz')
        summary = os.path.join(tmpdir, 'summary.csv')
        report = os.path.join(tmpdir, 'report.txt')
        contig_report = os.path.join(tmpdir, 'contigs.csv')

        args = build_args(genome=genome, matches_sig=matches_sig,
                          lineages_csv=lineages_csv, clean=clean, dirty=dirty,
                          summary=summary, report=report,
                          contig_report=contig_report, match_rank='order')

        main(args)
        
        assert os.path.exists(clean)
        assert os.path.exists(dirty)
        assert os.path.exists(report)
        assert os.path.exists(summary)
        assert os.path.exists(contig_report)

        with open(report, 'rt') as fp:
            data = fp.read()
            assert 'All genome hashes belong to one lineage! Nothing to do.' in data
