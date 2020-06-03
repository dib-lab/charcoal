"Test main() behavior, including file creation etc."

import types
import tempfile
import shutil
import os.path
import random


from charcoal.just_taxonomy import main
from .test_decontam import load_first_chunk


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


def build_args_with_tmpdir(tmpdir, **kwargs):
    args = types.SimpleNamespace(**kwargs)

    if 'clean' not in kwargs:
        args.clean = os.path.join(tmpdir, 'clean.fa.gz')
    if 'dirty' not in kwargs:
        args.dirty = os.path.join(tmpdir, 'dirty.fa.gz')
    if 'summary' not in kwargs:
        args.summary = os.path.join(tmpdir, 'summary.csv')
    if 'report' not in kwargs:
        args.report = os.path.join(tmpdir, 'report.txt')
    if 'contig_report' not in kwargs:
        args.contig_report = os.path.join(tmpdir, 'contigs.csv')

    return args


def test_main_1():
    # basic test - run it!
    genome = get_test_data('genomes/2.fa.gz')
    matches_sig = get_test_data('2.fa.gz.gather-matches.sig.gz')
    lineages_csv = get_test_data('test-match-lineages.csv')

    with TempDirectory() as tmpdir:
        args = build_args_with_tmpdir(tmpdir,
                                      genome=genome, matches_sig=matches_sig,
                                      lineages_csv=lineages_csv,
                                      match_rank='order')

        main(args)

        assert os.path.exists(args.clean)
        assert os.path.exists(args.dirty)
        assert os.path.exists(args.report)
        assert os.path.exists(args.summary)
        assert os.path.exists(args.contig_report)

        with open(args.report, 'rt') as fp:
            data = fp.read()
            assert 'All genome hashes belong to one lineage! Nothing to do.' in data


def test_main_2_unknown_lineage():
    # no matches
    genome = get_test_data('genomes/2.fa.gz')
    lineages_csv = get_test_data('test-match-lineages.csv')

    with TempDirectory() as tmpdir:
        matches_sig = os.path.join(tmpdir, 'empty.sig')
        # create empty file
        with open(matches_sig, 'wt') as fp:
            pass

        args = build_args_with_tmpdir(tmpdir,
                                      genome=genome, matches_sig=matches_sig,
                                      lineages_csv=lineages_csv,
                                      match_rank='order')

        main(args)

        assert os.path.exists(args.clean)
        assert os.path.exists(args.dirty)
        assert os.path.exists(args.report)
        assert os.path.exists(args.summary)
        assert os.path.exists(args.contig_report)

        with open(args.report, 'rt') as fp:
            data = fp.read()
            assert 'no matches to this genome were found in the database; nothing to do' in data


def test_main_3_unknown_lineage():
    # ???
    chunk = load_first_chunk(get_test_data('genomes/2.fa.gz'))[0]
    matches_sig = get_test_data('2.fa.gz.gather-matches.sig.gz')
    lineages_csv = get_test_data('test-match-lineages.csv')

    random_dna = ['A', 'C', 'G', 'T'] * int(50000 * 10)
    random.shuffle(random_dna)
    random_dna = "".join(random_dna)

    with TempDirectory() as tmpdir:
        genomefile = os.path.join(tmpdir, 'genome.fa')
        with open(genomefile, 'wt') as fp:
            print(f'>{chunk.name}\n{chunk.sequence}', file=fp)
            print(f'>random\n{random_dna}\n', file=fp)

        lineage = ''
        args = build_args_with_tmpdir(tmpdir,
                                      genome=genomefile,
                                      matches_sig=matches_sig,
                                      lineages_csv=lineages_csv,
                                      lineage=lineage,
                                      match_rank='order',
                                      force=False)

        main(args)

        assert os.path.exists(args.clean)
        assert os.path.exists(args.dirty)
        assert os.path.exists(args.report)
        assert os.path.exists(args.summary)
        assert os.path.exists(args.contig_report)

        with open(args.report, 'rt') as fp:
            data = fp.read()
            print(data)
            assert 'ERROR: fraction of total identified hashes (f_ident) < 10%.' in data


def test_main_4_unknown_lineage():
    # ???
    chunk = load_first_chunk(get_test_data('genomes/2.fa.gz'))[0]
    matches_sig = get_test_data('2.fa.gz.gather-matches.sig.gz')
    lineages_csv = get_test_data('test-match-lineages.csv')

    random_dna = ['A', 'C', 'G', 'T'] * int(50000 * 10)
    random.shuffle(random_dna)
    random_dna = "".join(random_dna)

    with TempDirectory() as tmpdir:
        genomefile = os.path.join(tmpdir, 'genome.fa')
        with open(genomefile, 'wt') as fp:
            print(f'>{chunk.name}\n{chunk.sequence}', file=fp)
            print(f'>random\n{random_dna}\n', file=fp)

        lineage = 'Bacteria;Verrucomicrobia;Verrucomicrobiae;Verrucomicrobiales'
        args = build_args_with_tmpdir(tmpdir,
                                      genome=genomefile,
                                      matches_sig=matches_sig,
                                      lineages_csv=lineages_csv,
                                      lineage=lineage,
                                      match_rank='order',
                                      force=False)

        main(args)

        assert os.path.exists(args.clean)
        assert os.path.exists(args.dirty)
        assert os.path.exists(args.report)
        assert os.path.exists(args.summary)
        assert os.path.exists(args.contig_report)

        with open(args.report, 'rt') as fp:
            data = fp.read()
            print(data)
            assert 'Using provided lineage as genome lineage.' in data
