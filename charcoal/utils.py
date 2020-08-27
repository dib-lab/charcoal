"""
utility functions for charcoal.
"""
import json
from collections import defaultdict, Counter

import sourmash
from sourmash.lca import lca_utils, LineagePair


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


def get_idents_for_hashval(lca_db, hashval):
    "Get the identifiers associated with this hashval."
    idx_list = lca_db.hashval_to_idx.get(hashval, [])
    for idx in idx_list:
        ident = lca_db.idx_to_ident[idx]
        yield ident


def gather_lca_assignments(hashvals, rank, dblist, ldb):
    """
    Collect lineage assignments from across all the databases for all the
    hashvals.
    """
    assignments = defaultdict(set)
    for hashval in hashvals:
        for lca_db in dblist:
            for ident in get_idents_for_hashval(lca_db, hashval):
                lineage = ldb.ident_to_lineage[ident]

                if rank:
                    lineage = pop_to_rank(lineage, rank)
                assignments[hashval].add(lineage)

    return assignments


def count_lca_for_assignments(assignments):
    """
    For each hashval, count the LCA across its assignments.
    """
    counts = Counter()
    for hashval in assignments:

        # for each list of tuple_info [(rank, name), ...] build
        # a tree that lets us discover lowest-common-ancestor.
        lineages = assignments[hashval]
        tree = sourmash.lca.build_tree(lineages)

        # now find either a leaf or the first node with multiple
        # children; that's our lowest-common-ancestor node.
        lca, reason = sourmash.lca.find_lca(tree)
        counts[lca] += 1

    return counts


def pretty_print_lineage(lin):
    "Nice output names for lineages."
    if not lin:
        return f'** no assignment **'
    elif lin[-1].rank == 'strain':
        strain = lin[-1].name
        return f'{strain}'
    elif lin[-1].rank == 'species':
        species = lin[-1].name
        return f'{species}'
    else:
        return f'{lin[-1].rank} {lin[-1].name}'


def pretty_print_lineage2(lin, rank):
    "Nice output names for lineages."
    if not lin:
        return f'** no assignment **'

    lin = pop_to_rank(lin, rank)
    return sourmash.lca.display_lineage(lin)


class WriteAndTrackFasta(object):
    def __init__(self, outfp, mh_ex):
        self.minhash = mh_ex.copy_and_clear()
        self.outfp = outfp
        self.n = 0
        self.bp = 0

    def write(self, record, no_write=False):
        if not no_write:
            self.outfp.write(f'>{record.name}\n{record.sequence}\n')
        self.minhash.add_sequence(record.sequence, force=True)
        self.n += 1
        self.bp += len(record.sequence)

    def close(self):
        self.outfp.close()


def gather_at_rank(mh, lca_db, lin_db, match_rank):
    "Run gather, and aggregate at given rank."
    import copy
    minhash = copy.copy(mh)
    query_sig = sourmash.SourmashSignature(minhash)

    # do the gather:
    counts = Counter()
    while 1:
        results = lca_db.gather(query_sig, threshold_bp=0)
        if not results:
            break

        (match, match_sig, _) = results[0]

        # retrieve lineage & pop to match_rank
        match_ident = get_ident(match_sig)
        match_lineage = lin_db.ident_to_lineage[match_ident]
        match_lineage = pop_to_rank(match_lineage, match_rank)

        # count at match_rank
        common = match_sig.minhash.count_common(query_sig.minhash)
        counts[match_lineage] += common

        # finish out gather algorithm!
        minhash.remove_many(match_sig.minhash.hashes)
        query_sig = sourmash.SourmashSignature(minhash)

    # return!
    for lin, count in counts.most_common():
        yield lin, count


def summarize_at_rank(lincounts, rank):
    newcounts = Counter()
    for lin, count in lincounts:
        lin = pop_to_rank(lin, rank)
        newcounts[lin] += count

    return newcounts.most_common()


def get_ident(sig):
    "Hack and slash identifiers."
    ident = sig.name()
    ident = ident.split()[0]
    ident = ident.split('.')[0]
    return ident

def load_contigs_gather_json(filename):
    # load contigs JSON file - @CTB
    with open(filename, 'rt') as fp:
        contigs_d = json.load(fp)
        for k in contigs_d:
            (size, v) = contigs_d[k]
            vv = []
            for (lin, count) in v:
                vv.append(([ LineagePair(*x) for x in lin ], count))
            contigs_d[k] = (size, vv)

    return contigs_d


def is_contig_contaminated(genome_lineage, contig_taxlist, rank, match_count_threshold):
    taxlist_at_rank = summarize_at_rank(contig_taxlist, rank)

    top_hit = None
    if contig_taxlist:
        top_hit, count = contig_taxlist[0]
        if count < match_count_threshold:
            top_hit = None

    is_bad = False
    if genome_lineage and top_hit and not is_lineage_match(genome_lineage, top_hit, rank):
        is_bad = True

        # rescue?
        for hit, count in contig_taxlist[1:]:
            if is_lineage_match(genome_lineage, hit, rank):
                is_bad = False

    return is_bad
