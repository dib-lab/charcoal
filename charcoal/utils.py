"""
utility functions for charcoal.
"""
import math
import numpy as np
from numpy import genfromtxt
import screed

from sourmash.lca import lca_utils


def load_matrix_csv(filename):
    mat = genfromtxt(filename, delimiter=',')
    return mat

def make_distance_matrix(mat, delete_empty=False):
    """
    Construct distance matrix from metagenome x hash matrices.
    """
    n_hashes = mat.shape[1]
    n_orig_hashes = n_hashes

    # go through and normalize all the sample-presence vectors for each hash;
    # track those with all 0s for later removal.
    to_delete = []
    for i in range(n_hashes):
        if sum(mat[:, i]):
            mat[:, i] /= math.sqrt(np.dot(mat[:, i], mat[:, i]))
        else:
            to_delete.append(i)

    if delete_empty:
        # remove all columns with zeros
        print('removing {} null presence vectors'.format(len(to_delete)))
        for row_n in reversed(to_delete):
            mat = np.delete(mat, row_n, 1)

        assert mat.shape[1] == n_hashes - len(to_delete)

    n_hashes = mat.shape[1]

    # construct distance matrix using angular distance
    D = np.zeros((n_hashes, n_hashes))
    for i in range(n_hashes):
        for j in range(n_hashes):
            cos_sim = np.dot(mat[:, i], mat[:, j])
            cos_sim = min(cos_sim, 1.0)
            ang_sim = 1 - 2*math.acos(cos_sim) / math.pi
            D[i][j] = ang_sim

    # done!
    return D, n_orig_hashes


def is_lineage_match(lin_a, lin_b, rank):
    """
    check to see if two lineages are a match down to given rank.
    """
    for a, b in zip(lin_a, lin_b):
        assert a.rank == b.rank
        if a.rank == rank:
            if a == b:
                return 1
        if a != b:
            return 0

    return 0


def pop_to_rank(lin, rank):
    "Remove lineage tuples from given lineage `lin` until `rank` is reached."
    lin = list(lin)

    txl = lca_utils.taxlist()
    before_rank = []
    for txl_rank in txl:
        if txl_rank != rank:
            before_rank.append(txl_rank)
        else:
            break

    # are we already above rank?
    if lin and lin[-1].rank in before_rank:
        return tuple(lin)

    while lin and lin[-1].rank != rank:
        lin.pop()

    return tuple(lin)


class HashesToTaxonomy(object):
    def __init__(self, genome_file, ksize, scaled, fragment_size, lca_db_file):
        self.genome_file = genome_file
        self.ksize = ksize
        self.scaled = scaled
        self.fragment_size = fragment_size
        self.lca_db_file = lca_db_file

        self.d = {}

    def __setitem__(self, hashval, lineage):
        self.d[hashval] = lineage

    def __getitem__(self, hashval):
        return self.d[hashval]

    def __len__(self):
        return len(self.d)

    def __iter__(self):
        return iter(self.d)

    def items(self):
        return self.d.items()


class HashesToLengths(object):
    def __init__(self, genome_file, ksize, scaled, fragment_size):
        self.genome_file = genome_file
        self.ksize = ksize
        self.scaled = scaled
        self.fragment_size = fragment_size

        self.d = {}

    def __setitem__(self, hashval, length):
        self.d[hashval] = length

    def __getitem__(self, hashval):
        return self.d[hashval]

    def __len__(self):
        return len(self.d)

    def __iter__(self):
        return iter(self.d)

    def items(self):
        return self.d.items()


class MetagenomesMatrix(object):
    def __init__(self, genome_file, query_hashlist, query_fragment_size, ksize):
        self.genome_file = genome_file
        self.query_hashlist = list(sorted(query_hashlist))
        self.query_fragment_size = query_fragment_size
        self.ksize = ksize
        self.mat = None


class GenomeShredder(object):
    def __init__(self, genome_file, fragment_size):
        self.genome_file = genome_file
        self.fragment_size = fragment_size

    def __iter__(self):
        fragment_size = self.fragment_size

        for record in screed.open(self.genome_file):
            if not fragment_size:
                yield record.name, record.sequence, 0, len(record.sequence)
            else:
                for start in range(0, len(record.sequence), fragment_size):
                    seq = record.sequence[start:start + fragment_size]
                    yield record.name, seq, start, start + len(seq)
