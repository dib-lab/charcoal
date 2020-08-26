"""Lineage database and associated utilities.

Extracted from sourmash LCA Databases.

TODO:
* implement load, save

CTB:
* add strict taxonomy enforcement within? ldb-specific taxlist?
"""

from __future__ import print_function, division
import json
import gzip
from collections import OrderedDict, defaultdict
import functools
import pytest

from sourmash.lca import LineagePair


def cached_property(fun):
    """A memoize decorator for class properties."""
    @functools.wraps(fun)
    def get(self):
        try:
            return self._cache[fun]
        except AttributeError:
            self._cache = {}
        except KeyError:
            pass
        ret = self._cache[fun] = fun(self)
        return ret
    return property(get)


class LineageDB(object):
    """
    An in-memory database for taxonomic lineages.

    Identifiers `ident` must be unique.

    Integer `lid` indices can be used as keys in dictionary attributes:
    * `lid_to_lineage`, to get a lineage for that `lid`.
    * `lid_to_idents`, to get set of identifiers for that `lid`.

    `lineage_to_lid` is a dictionary with tuples of LineagePair as keys,
    `lid` as values.

    `ident_to_lid` is a dictionary from unique str identifer to integer `lid`.
    `ident_to_lineage` is a dict from identifier to lineage.
    """
    def __init__(self):
        self.lineage_to_idents = defaultdict(set)
        self.ident_to_lineage = {}

    def _invalidate_cache(self):
        if hasattr(self, '_cache'):
            del self._cache

    def insert(self, ident, lineage):
        """Add a new identity / lineage pair into the database.

        'ident' must be a unique string identifer across this database.

        'lineage', if specified, must contain a tuple of LineagePair objects.
        """
        if ident in self.ident_to_lineage:
            raise ValueError("identifier {} is already in this lineage db.".format(ident))

        # before adding, invalide any caching from @cached_property
        self._invalidate_cache()

        try:
            lineage = tuple(lineage)
            self.lineage_to_idents[lineage].add(ident)
        except TypeError:
            raise ValueError('lineage cannot be used as a key?!')

        self.ident_to_lineage[ident] = lineage

    def __repr__(self):
        return "LineageDB('{}')".format(self.filename)

    @classmethod
    def load(cls, db_name):
        "Load LCA_Database from a JSON file."
        from .lca_utils import taxlist, LineagePair

        xopen = open
        if db_name.endswith('.gz'):
            xopen = gzip.open

        with xopen(db_name, 'rt') as fp:
            load_d = {}
            try:
                load_d = json.load(fp)
            except json.decoder.JSONDecodeError:
                pass

            if not load_d:
                raise ValueError("cannot parse database file '{}' as JSON; invalid format.")

            version = None
            db_type = None
            try:
                version = load_d.get('version')
                db_type = load_d.get('type')
            except AttributeError:
                pass

            if db_type != 'sourmash_lca':
                raise ValueError("database file '{}' is not an LCA db.".format(db_name))

            if version != '2.0' or 'lid_to_lineage' not in load_d:
                raise ValueError("Error! This is an old-style LCA DB. You'll need to rebuild or download a newer one.")

            ksize = int(load_d['ksize'])
            scaled = int(load_d['scaled'])

            db = cls(ksize, scaled)

            # convert lineage_dict to proper lineages (tuples of LineagePairs)
            lid_to_lineage_2 = load_d['lid_to_lineage']
            lid_to_lineage = {}
            lineage_to_lid = {}
            for k, v in lid_to_lineage_2.items():
                v = dict(v)
                vv = []
                for rank in taxlist():
                    name = v.get(rank, '')
                    vv.append(LineagePair(rank, name))

                vv = tuple(vv)
                lid_to_lineage[int(k)] = vv
                lineage_to_lid[vv] = int(k)
            db.lid_to_lineage = lid_to_lineage
            db.lineage_to_lid = lineage_to_lid

            # convert hashval -> lineage index keys to integers (looks like
            # JSON doesn't have a 64 bit type so stores them as strings)
            hashval_to_idx_2 = load_d['hashval_to_idx']
            hashval_to_idx = {}

            for k, v in hashval_to_idx_2.items():
                hashval_to_idx[int(k)] = v
            db.hashval_to_idx = hashval_to_idx

            db.ident_to_name = load_d['ident_to_name']
            db.ident_to_idx = load_d['ident_to_idx']

            db.idx_to_lid = {}
            for k, v in load_d['idx_to_lid'].items():
                db.idx_to_lid[int(k)] = v

        db.filename = db_name

        return db

    def save(self, db_name):
        "Save LCA_Database to a JSON file."
        xopen = open
        if db_name.endswith('.gz'):
            xopen = gzip.open

        with xopen(db_name, 'wt') as fp:
            # use an OrderedDict to preserve output order
            save_d = OrderedDict()
            save_d['version'] = '2.0'
            save_d['type'] = 'sourmash_lca'
            save_d['license'] = 'CC0'
            save_d['ksize'] = self.ksize
            save_d['scaled'] = self.scaled

            # convert lineage internals from tuples to dictionaries
            d = OrderedDict()
            for k, v in self.lid_to_lineage.items():
                d[k] = dict([ (vv.rank, vv.name) for vv in v ])
            save_d['lid_to_lineage'] = d

            # convert values from sets to lists, so that JSON knows how to save
            save_d['hashval_to_idx'] = \
               dict((k, list(v)) for (k, v) in self.hashval_to_idx.items())

            save_d['ident_to_name'] = self.ident_to_name
            save_d['ident_to_idx'] = self.ident_to_idx
            save_d['idx_to_lid'] = self.idx_to_lid
            save_d['lid_to_lineage'] = self.lid_to_lineage
            
            json.dump(save_d, fp)


def test_lineage_db_1():
    lineage = ((LineagePair('rank1', 'name1'),
                LineagePair('rank2', 'name2')))

    ldb = LineageDB()
    ldb.insert('uniq', lineage)

    assert 'uniq' in ldb.lineage_to_idents[lineage]
    assert ldb.ident_to_lineage['uniq'] == lineage


def test_lineage_db_1_tuple():
    # list here cannot be tuple-ized for use as a lineage
    lineage = ([LineagePair('rank1', 'name1'),
                LineagePair('rank2', 'name2')],)

    ldb = LineageDB()

    with pytest.raises(ValueError):
        ldb.insert('uniq', lineage)


def test_lineage_db_1_non_iter():
    # try a non-iterable => fail.
    lineage = 1

    ldb = LineageDB()

    with pytest.raises(ValueError):
        ldb.insert('uniq', lineage)
