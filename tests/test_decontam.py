from charcoal import just_taxonomy, utils
from sourmash.lca import LineagePair, taxlist


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
