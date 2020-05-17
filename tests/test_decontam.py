from io import StringIO
import gzip

import screed

from charcoal import just_taxonomy, utils
from charcoal.lineage_db import LineageDB

import sourmash
from sourmash.lca import LCA_Database
from sourmash.lca import LineagePair, taxlist
from sourmash.lca.command_index import load_taxonomy_assignments


def make_lineage(lineage):
    "Turn a ; or ,-separated set of lineages into a tuple of LineagePair objs."
    lin = lineage.split(';')
    if len(lin) == 1:
        lin = lineage.split(',')
    lin = [ LineagePair(rank, n) for (rank, n) in zip(taxlist(), lin) ]

    return lin


class ReportingCapture(object):
    "Capture reporting output via 'report' fn."
    def __init__(self):
        self.lines = []

    def __call__(self, s):
        self.lines.append(s)

    def contains(self, s_try):
        for x in self.lines:
            if s_try in x:
                return True
        return False


class FakeFastaWriter(object):
    "Test fixture to record FASTA output from ContigDecontaminator class."
    def __init__(self):
        self.n = 0
        self.bp = 0
        self.names = []

    def write(self, record):
        self.names.append(record.name)
        self.bp += len(record.sequence)
        self.n += 1


def make_lca_and_lineages(match_files, lineages_csv, scaled, ksize,
                          start_column=3):
    "Build an LCA_Database and LineageDB."
    lca_db = LCA_Database(ksize=ksize, scaled=scaled)
    lin_db = LineageDB()

    # load matches into a list
    siglist = []
    for filename in match_files:
        siglist += list(sourmash.load_signatures(filename))

    # pull in taxonomic assignments
    tax_assign, _ = load_taxonomy_assignments(lineages_csv,
                                              start_column=start_column)

    # build database of matches & lineages!
    for ss in siglist:
        print(ss.name(), ss.minhash.scaled)
        ident = just_taxonomy.get_ident(ss)
        lineage = tax_assign[ident]

        lca_db.insert(ss, ident=ident)
        lin_db.insert(ident, lineage)

    return lca_db, lin_db


def load_first_chunk(filename, chunksize=50000):
    # load a first small chunk from a contiguous genome. eventually,
    # need to fix to pay attention to contigs...
    with gzip.open(filename, 'rb') as fp:
        seqname = fp.readline().strip()[1:]
        seqdata = fp.read(chunksize).decode('utf-8')
        seqdata = seqdata.split()
        seqdata = "".join(seqdata)

    return [screed.Record(seqname, seqdata)]


def test_choose_genome_lineage_1():
    # unambiguous choice - good lca lineage, nothing else provided

    lca_lin = "Bacteria;Proteobacteria;Gammaproteobacteria;Alteromonadales;Shewanellaceae;Shewanella;Shewanella baltica;Shewanella baltica OS185"
    lca_lin = make_lineage(lca_lin)
    provided_lin = ""
    f_ident = 1.0
    f_major = 1.0
    x = just_taxonomy.choose_genome_lineage(lca_lin, provided_lin,
                                            f_ident, f_major,
                                            ReportingCapture())
    (chosen_lin, comment) = x

    assert chosen_lin == utils.pop_to_rank(lca_lin, 'genus')
    assert not comment


def test_choose_genome_lineage_2():
    # unambiguous choice - good provided lineage

    lca_lin = "Bacteria;Proteobacteria;Gammaproteobacteria;Alteromonadales;Shewanellaceae;Shewanella;Shewanella baltica"
    lca_lin = make_lineage(lca_lin)

    provided_lin = "Bacteria;Verrucomicrobia;Verrucomicrobiae;Verrucomicrobiales;Akkermansiaceae;Akkermansia;Akkermansia muciniphila"
    provided_lin = make_lineage(provided_lin)
    assert lca_lin != provided_lin

    f_ident = 1.0
    f_major = 1.0
    x = just_taxonomy.choose_genome_lineage(lca_lin, provided_lin,
                                            f_ident, f_major,
                                            ReportingCapture())
    (chosen_lin, comment) = x

    assert chosen_lin == utils.pop_to_rank(provided_lin, 'genus')
    assert not comment


def test_choose_genome_lineage_3():
    # agreement between lca and provided

    lca_lin = "Bacteria;Proteobacteria;Gammaproteobacteria;Alteromonadales;Shewanellaceae;Shewanella;Shewanella baltica"
    lca_lin = make_lineage(lca_lin)

    provided_lin = "Bacteria;Proteobacteria;Gammaproteobacteria;Alteromonadales;Shewanellaceae;Shewanella;Shewanella baltica;Shewanella baltica OS185"
    provided_lin = make_lineage(provided_lin)

    f_ident = 1.0
    f_major = 1.0
    reporting = ReportingCapture()
    x = just_taxonomy.choose_genome_lineage(lca_lin, provided_lin,
                                            f_ident, f_major, reporting)
    chosen_lin, comment = x

    assert reporting.contains('provided lineage agrees with k-mer class')

    assert chosen_lin == utils.pop_to_rank(lca_lin, 'genus')
    assert not comment


def test_choose_genome_lineage_4():
    # disagreement between lca and provided

    lca_lin = "Bacteria;Proteobacteria;Gammaproteobacteria;Alteromonadales;Shewanellaceae;Shewanella;Shewanella baltica"
    lca_lin = make_lineage(lca_lin)

    provided_lin = "Bacteria;Verrucomicrobia;Verrucomicrobiae;Verrucomicrobiales;Akkermansiaceae;Akkermansia;Akkermansia muciniphila"
    provided_lin = make_lineage(provided_lin)

    f_ident = 1.0
    f_major = 1.0
    reporting = ReportingCapture()
    x = just_taxonomy.choose_genome_lineage(lca_lin, provided_lin,
                                            f_ident, f_major, reporting)
    chosen_lin, comment = x

    assert reporting.contains('provided lineage disagrees with k-mer class')

    assert chosen_lin == utils.pop_to_rank(provided_lin, 'genus')
    assert not comment


def test_choose_genome_lineage_5():
    # bad choice - bad lca lineage, nothing else provided

    lca_lin = "Bacteria;Proteobacteria;Gammaproteobacteria;Alteromonadales;Shewanellaceae;Shewanella;Shewanella baltica;Shewanella baltica OS185"
    lca_lin = make_lineage(lca_lin)
    provided_lin = ""
    f_ident = 0.05
    f_major = 1.0
    reporting = ReportingCapture()
    x = just_taxonomy.choose_genome_lineage(lca_lin, provided_lin,
                                            f_ident, f_major, reporting)
    (chosen_lin, comment) = x

    print(reporting.lines)
    assert reporting.contains('Please provide a lineage for this genome')
    assert 'too few identifiable hashes' in comment


def test_choose_genome_lineage_6():
    # bad choice - bad lca lineage, nothing else provided

    lca_lin = "Bacteria;Proteobacteria;Gammaproteobacteria;Alteromonadales;Shewanellaceae;Shewanella;Shewanella baltica;Shewanella baltica OS185"
    lca_lin = make_lineage(lca_lin)
    provided_lin = ""
    f_ident = 1.0
    f_major = 0.15
    reporting = ReportingCapture()
    x = just_taxonomy.choose_genome_lineage(lca_lin, provided_lin,
                                            f_ident, f_major, reporting)
    (chosen_lin, comment) = x

    print(reporting.lines)
    assert reporting.contains('Please provide a lineage for this genome')
    assert 'too few hashes in major lineage' in comment


def test_cleaner_class_1():
    # your basic clean genome test
    genome_lineage = "Bacteria;Proteobacteria;Gammaproteobacteria;Alteromonadales;Shewanellaceae;Shewanella;Shewanella baltica;Shewanella baltica OS185"
    genome_lineage = make_lineage(genome_lineage)

    match_rank = 'genus'
    empty_mh = sourmash.MinHash(n=0, ksize=31, scaled=10000)

    matches_file = 'tests/test-data/63.fa.gz.gather-matches.sig.gz'
    lineages_csv = 'test-data/test-match-lineages.csv'
    lca_db, lin_db = make_lca_and_lineages([matches_file], lineages_csv,
                                           empty_mh.scaled, empty_mh.ksize)

    cleaner = just_taxonomy.ContigsDecontaminator(genome_lineage,
                                                  match_rank,
                                                  empty_mh,
                                                  lca_db, lin_db)
    cleaner.clean_out = FakeFastaWriter()
    cleaner.dirty_out = FakeFastaWriter()

    inp_iter = load_first_chunk('test-data/genomes/63.fa.gz')

    report_fp = StringIO()
    cleaner.clean_contigs(inp_iter, report_fp)

    assert cleaner.dirty_out.n == 0
    assert cleaner.clean_out.n == 1
    assert cleaner.clean_out.bp == 49383
    assert cleaner.clean_out.names == [b'NC_011663.1 Shewanella baltica OS223, complete genome']


def test_cleaner_class_2():
    # specify wrong lineage - everything must go :)
    genome_lineage = "Bacteria;Verrucomicrobia;Verrucomicrobiae;Verrucomicrobiales;Akkermansiaceae;Akkermansia;Akkermansia muciniphila"
    genome_lineage = make_lineage(genome_lineage)

    match_rank = 'genus'
    empty_mh = sourmash.MinHash(n=0, ksize=31, scaled=10000)

    matches_file = 'tests/test-data/63.fa.gz.gather-matches.sig.gz'
    lineages_csv = 'test-data/test-match-lineages.csv'
    lca_db, lin_db = make_lca_and_lineages([matches_file], lineages_csv,
                                           empty_mh.scaled, empty_mh.ksize)

    cleaner = just_taxonomy.ContigsDecontaminator(genome_lineage,
                                                  match_rank,
                                                  empty_mh,
                                                  lca_db, lin_db)
    cleaner.clean_out = FakeFastaWriter()
    cleaner.dirty_out = FakeFastaWriter()

    inp_iter = load_first_chunk('test-data/genomes/63.fa.gz')

    report_fp = StringIO()
    cleaner.clean_contigs(inp_iter, report_fp)

    assert cleaner.clean_out.n == 0
    assert cleaner.dirty_out.n == 1
    assert cleaner.dirty_out.bp == 49383
    assert cleaner.dirty_out.names == [b'NC_011663.1 Shewanella baltica OS223, complete genome']


def test_cleaner_class_3():
    # specify two contigs, with mixed lineage - one clean one dirty
    genome_lineage = "Bacteria;Verrucomicrobia;Verrucomicrobiae;Verrucomicrobiales;Akkermansiaceae;Akkermansia;Akkermansia muciniphila"
    genome_lineage = make_lineage(genome_lineage)

    match_rank = 'genus'
    empty_mh = sourmash.MinHash(n=0, ksize=31, scaled=10000)

    matches_file1 = 'tests/test-data/2.fa.gz.gather-matches.sig.gz'
    matches_file2 = 'tests/test-data/63.fa.gz.gather-matches.sig.gz'
    lineages_csv = 'test-data/test-match-lineages.csv'
    lca_db, lin_db = make_lca_and_lineages([matches_file1, matches_file2],
                                           lineages_csv,
                                           empty_mh.scaled, empty_mh.ksize)

    cleaner = just_taxonomy.ContigsDecontaminator(genome_lineage,
                                                  match_rank,
                                                  empty_mh,
                                                  lca_db, lin_db)
    cleaner.clean_out = FakeFastaWriter()
    cleaner.dirty_out = FakeFastaWriter()

    chunk1 = load_first_chunk('test-data/genomes/63.fa.gz')
    chunk2 = load_first_chunk('test-data/genomes/2.fa.gz')
    inp_iter = chunk1 + chunk2

    report_fp = StringIO()
    cleaner.clean_contigs(inp_iter, report_fp)

    assert cleaner.dirty_out.n == 1
    assert cleaner.dirty_out.bp == 49383
    assert cleaner.dirty_out.names == [b'NC_011663.1 Shewanella baltica OS223, complete genome']

    assert cleaner.clean_out.n == 1
    assert cleaner.clean_out.bp == 49296
    assert cleaner.clean_out.names == [b'CP001071.1 Akkermansia muciniphila ATCC BAA-835, complete genome']

    # double-check lineage juuuuust to confirm
    clean_name = cleaner.clean_out.names[0]
    clean_name = clean_name.decode('utf-8')
    ident = clean_name.split('.')[0]
    clean_lineage = lin_db.ident_to_lineage[ident]
    assert utils.is_lineage_match(clean_lineage, genome_lineage, 'genus')


def test_cleaner_class_4():
    # specify two bacterial contigs, and ask for superkingdom level cleaning.
    # (should be nothing dirty)
    genome_lineage = "Bacteria;Verrucomicrobia;Verrucomicrobiae;Verrucomicrobiales;Akkermansiaceae;Akkermansia;Akkermansia muciniphila"
    genome_lineage = make_lineage(genome_lineage)

    match_rank = 'superkingdom'
    empty_mh = sourmash.MinHash(n=0, ksize=31, scaled=10000)

    matches_file1 = 'tests/test-data/2.fa.gz.gather-matches.sig.gz'
    matches_file2 = 'tests/test-data/63.fa.gz.gather-matches.sig.gz'
    lineages_csv = 'test-data/test-match-lineages.csv'
    lca_db, lin_db = make_lca_and_lineages([matches_file1, matches_file2],
                                           lineages_csv,
                                           empty_mh.scaled, empty_mh.ksize)

    cleaner = just_taxonomy.ContigsDecontaminator(genome_lineage,
                                                  match_rank,
                                                  empty_mh,
                                                  lca_db, lin_db)
    cleaner.clean_out = FakeFastaWriter()
    cleaner.dirty_out = FakeFastaWriter()

    chunk1 = load_first_chunk('test-data/genomes/63.fa.gz')
    chunk2 = load_first_chunk('test-data/genomes/2.fa.gz')
    inp_iter = chunk1 + chunk2

    report_fp = StringIO()
    cleaner.clean_contigs(inp_iter, report_fp)

    assert cleaner.dirty_out.n == 0
    assert cleaner.dirty_out.bp == 0
    assert cleaner.dirty_out.names == []

    assert cleaner.clean_out.n == 2
    assert cleaner.clean_out.bp == 49296 + 49383
    assert cleaner.clean_out.names == [b'NC_011663.1 Shewanella baltica OS223, complete genome', b'CP001071.1 Akkermansia muciniphila ATCC BAA-835, complete genome']

    # double-check lineage juuuuust to confirm
    for clean_name in cleaner.clean_out.names:
        clean_name = clean_name.decode('utf-8')
        ident = clean_name.split('.')[0]
        clean_lineage = lin_db.ident_to_lineage[ident]
        assert utils.is_lineage_match(clean_lineage, genome_lineage,
                                      match_rank)


def test_cleaner_class_5():
    # specify two bacterial contigs, and ask for phylum level cleaning.
    # => one clean one dirty
    genome_lineage = "Bacteria;Verrucomicrobia;Verrucomicrobiae;Verrucomicrobiales;Akkermansiaceae;Akkermansia;Akkermansia muciniphila"
    genome_lineage = make_lineage(genome_lineage)

    match_rank = 'genus'
    empty_mh = sourmash.MinHash(n=0, ksize=31, scaled=10000)

    matches_file1 = 'tests/test-data/2.fa.gz.gather-matches.sig.gz'
    matches_file2 = 'tests/test-data/63.fa.gz.gather-matches.sig.gz'
    lineages_csv = 'test-data/test-match-lineages.csv'
    lca_db, lin_db = make_lca_and_lineages([matches_file1, matches_file2],
                                           lineages_csv,
                                           empty_mh.scaled, empty_mh.ksize)

    cleaner = just_taxonomy.ContigsDecontaminator(genome_lineage,
                                                  match_rank,
                                                  empty_mh,
                                                  lca_db, lin_db)
    cleaner.clean_out = FakeFastaWriter()
    cleaner.dirty_out = FakeFastaWriter()

    chunk1 = load_first_chunk('test-data/genomes/63.fa.gz')
    chunk2 = load_first_chunk('test-data/genomes/2.fa.gz')
    inp_iter = chunk1 + chunk2

    report_fp = StringIO()
    cleaner.clean_contigs(inp_iter, report_fp)

    assert cleaner.dirty_out.n == 1
    assert cleaner.dirty_out.bp == 49383
    assert cleaner.dirty_out.names == [b'NC_011663.1 Shewanella baltica OS223, complete genome']

    assert cleaner.clean_out.n == 1
    assert cleaner.clean_out.bp == 49296
    assert cleaner.clean_out.names == [b'CP001071.1 Akkermansia muciniphila ATCC BAA-835, complete genome']

    # double-check lineage juuuuust to confirm
    clean_name = cleaner.clean_out.names[0]
    clean_name = clean_name.decode('utf-8')
    ident = clean_name.split('.')[0]
    clean_lineage = lin_db.ident_to_lineage[ident]
    assert utils.is_lineage_match(clean_lineage, genome_lineage, 'genus')


def test_cleaner_gather_method_1():
    # specify two bacterial contigs, and ask for phylum level cleaning.
    # => one clean one dirty
    genome_lineage = "Bacteria;Verrucomicrobia;Verrucomicrobiae;Verrucomicrobiales;Akkermansiaceae;Akkermansia;Akkermansia muciniphila"
    genome_lineage = make_lineage(genome_lineage)

    match_rank = 'genus'
    empty_mh = sourmash.MinHash(n=0, ksize=31, scaled=10000)

    matches_file1 = 'tests/test-data/2.fa.gz.gather-matches.sig.gz'
    matches_file2 = 'tests/test-data/63.fa.gz.gather-matches.sig.gz'
    lineages_csv = 'test-data/test-match-lineages.csv'
    lca_db, lin_db = make_lca_and_lineages([matches_file1, matches_file2],
                                           lineages_csv,
                                           empty_mh.scaled, empty_mh.ksize)

    cleaner = just_taxonomy.ContigsDecontaminator(genome_lineage,
                                                  match_rank,
                                                  empty_mh,
                                                  lca_db, lin_db)
    report_fp = StringIO()

    chunk1 = load_first_chunk('test-data/genomes/63.fa.gz')

    mh = empty_mh.copy_and_clear()
    mh.add_sequence(chunk1[0].sequence, force=True)
    is_clean = cleaner.check_gather(chunk1[0], mh, report_fp)
    assert not is_clean

    chunk2 = load_first_chunk('test-data/genomes/2.fa.gz')

    mh = empty_mh.copy_and_clear()
    mh.add_sequence(chunk2[0].sequence, force=True)
    is_clean = cleaner.check_gather(chunk2[0], mh, report_fp)
    assert is_clean


def test_cleaner_gather_method_1():
    # specify two bacterial contigs, and ask for phylum level cleaning.
    # => one clean one dirty
    genome_lineage = "Bacteria;Verrucomicrobia;Verrucomicrobiae;Verrucomicrobiales;Akkermansiaceae;Akkermansia;Akkermansia muciniphila"
    genome_lineage = make_lineage(genome_lineage)

    match_rank = 'genus'
    empty_mh = sourmash.MinHash(n=0, ksize=31, scaled=10000)

    matches_file1 = 'tests/test-data/2.fa.gz.gather-matches.sig.gz'
    matches_file2 = 'tests/test-data/63.fa.gz.gather-matches.sig.gz'
    lineages_csv = 'test-data/test-match-lineages.csv'
    lca_db, lin_db = make_lca_and_lineages([matches_file1, matches_file2],
                                           lineages_csv,
                                           empty_mh.scaled, empty_mh.ksize)

    cleaner = just_taxonomy.ContigsDecontaminator(genome_lineage,
                                                  match_rank,
                                                  empty_mh,
                                                  lca_db, lin_db)
    report_fp = StringIO()

    chunk1 = load_first_chunk('test-data/genomes/63.fa.gz')

    mh = empty_mh.copy_and_clear()
    mh.add_sequence(chunk1[0].sequence, force=True)
    is_clean = cleaner.check_gather(chunk1[0], mh, report_fp)
    assert not is_clean
    assert cleaner.n_reason_3 == 1

    chunk2 = load_first_chunk('test-data/genomes/2.fa.gz')

    mh = empty_mh.copy_and_clear()
    mh.add_sequence(chunk2[0].sequence, force=True)
    is_clean = cleaner.check_gather(chunk2[0], mh, report_fp)
    assert is_clean
    assert cleaner.n_reason_3 == 1        # should not be incremented


def test_cleaner_lca_method2_1():
    # specify two bacterial contigs, and ask for genus level cleaning.
    # => one clean one dirty, for reason 2
    genome_lineage = "Bacteria;Verrucomicrobia;Verrucomicrobiae;Verrucomicrobiales;Akkermansiaceae;Akkermansia;Akkermansia muciniphila"
    genome_lineage = make_lineage(genome_lineage)

    match_rank = 'genus'
    empty_mh = sourmash.MinHash(n=0, ksize=31, scaled=10000)

    matches_file1 = 'tests/test-data/2.fa.gz.gather-matches.sig.gz'
    matches_file2 = 'tests/test-data/63.fa.gz.gather-matches.sig.gz'
    lineages_csv = 'test-data/test-match-lineages.csv'
    lca_db, lin_db = make_lca_and_lineages([matches_file1, matches_file2],
                                           lineages_csv,
                                           empty_mh.scaled, empty_mh.ksize)

    cleaner = just_taxonomy.ContigsDecontaminator(genome_lineage,
                                                  match_rank,
                                                  empty_mh,
                                                  lca_db, lin_db)
    report_fp = StringIO()

    chunk1 = load_first_chunk('test-data/genomes/63.fa.gz')

    mh = empty_mh.copy_and_clear()
    mh.add_sequence(chunk1[0].sequence, force=True)
    is_clean = cleaner.check_lca(chunk1[0], mh, report_fp)
    assert not is_clean
    assert cleaner.n_reason_2 == 1

    chunk2 = load_first_chunk('test-data/genomes/2.fa.gz')

    mh = empty_mh.copy_and_clear()
    mh.add_sequence(chunk2[0].sequence, force=True)
    is_clean = cleaner.check_lca(chunk2[0], mh, report_fp)
    assert is_clean
    assert cleaner.n_reason_2 == 1        # should not be incremented


def test_cleaner_lca_method1_1():
    # label hashes in one contig with two different bacterial lineages
    # => one dirty, for reason 1
    genome_lineage = "Bacteria;Verrucomicrobia;Verrucomicrobiae;Verrucomicrobiales;Akkermansiaceae;Akkermansia;Akkermansia muciniphila"
    genome_lineage = make_lineage(genome_lineage)

    match_rank = 'genus'
    empty_mh = sourmash.MinHash(n=0, ksize=31, scaled=10000)

    matches_file1 = 'tests/test-data/2.fa.gz.gather-matches.sig.gz'
    lineages_csv = 'test-data/test-match-lineages.csv'
    lca_db, lin_db = make_lca_and_lineages([matches_file1],
                                           lineages_csv,
                                           empty_mh.scaled, empty_mh.ksize)

    # now, hack the lineages so that it looks like the same sequence has
    # a completely different lineage...
    matches_2 = sourmash.load_signatures('tests/test-data/2.fa.gz.gather-matches.sig.gz')
    matches_2 = list(matches_2)
    fake_lin = 'Bacteria,Proteobacteria,Gammaproteobacteria,Alteromonadales,Shewanellaceae,Shewanella,Shewanella baltica'
    fake_lin = make_lineage(fake_lin)

    lin_db.insert('NC_011663', fake_lin)
    for ss in matches_2:
        lca_db.insert(ss, 'NC_011663')

    # ok! now, go for cleaning...
    cleaner = just_taxonomy.ContigsDecontaminator(genome_lineage,
                                                  match_rank,
                                                  empty_mh,
                                                  lca_db, lin_db)
    report_fp = StringIO()

    chunk1 = load_first_chunk('test-data/genomes/2.fa.gz')

    mh = empty_mh.copy_and_clear()
    mh.add_sequence(chunk1[0].sequence, force=True)
    is_clean = cleaner.check_lca(chunk1[0], mh, report_fp,
                                 force_report=True)
    assert not is_clean
    assert cleaner.n_reason_1 == 1
