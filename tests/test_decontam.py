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
   lin = lineage.split(';')
   lin = [ LineagePair(rank, n) for (rank, n) in zip(taxlist(), lin) ]

   return lin


class ReportingCapture(object):
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
    def __init__(self):
        self.n = 0
        self.bp = 0
        self.names = []

    def write(self, record):
        self.names.append(record.name)
        self.bp += len(record.sequence)
        self.n += 1


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
    genome_lineage = "Bacteria;Proteobacteria;Gammaproteobacteria;Alteromonadales;Shewanellaceae;Shewanella;Shewanella baltica;Shewanella baltica OS185"
    genome_lineage = make_lineage(genome_lineage)

    match_rank = 'genus'
    empty_mh = sourmash.MinHash(n=0, ksize=31, scaled=10000)
    lca_db = LCA_Database(ksize=empty_mh.ksize, scaled=empty_mh.scaled)
    lin_db = LineageDB()

    siglist = list(sourmash.load_signatures('tests/test-data/63.fa.gz.gather-matches.sig'))
    tax_assign, _ = load_taxonomy_assignments('test-data/test-match-lineages.csv',
                                              start_column=3)

    # CTB: move to util or convenience function
    for ss in siglist:
        print(ss.name(), ss.minhash.scaled)
        ident = just_taxonomy.get_ident(ss)
        lineage = tax_assign[ident]

        lca_db.insert(ss, ident=ident)
        lin_db.insert(ident, lineage)


    cleaner = just_taxonomy.ContigsDecontaminator(genome_lineage,
                                                  match_rank,
                                                  empty_mh,
                                                  lca_db, lin_db)
    cleaner.clean_out = FakeFastaWriter()
    cleaner.dirty_out = FakeFastaWriter()

    # load a 50kb chunk of 63.fa.gz
    with gzip.open('test-data/genomes/63.fa.gz', 'rb') as fp:
        seqname = fp.readline().strip()
        seqdata = fp.read(50000).decode('utf-8')
        seqdata = seqdata.split()
        seqdata = "".join(seqdata)

    x = [screed.Record(seqname, seqdata)]

    report_fp = StringIO()
    cleaner.clean_contigs(x, report_fp)

    assert cleaner.clean_out.n == 1
    assert cleaner.clean_out.bp == 49383
    assert cleaner.clean_out.names == [b'>NC_011663.1 Shewanella baltica OS223, complete genome']
