__author__ = "alvaro barbeira"

import pandas
import copy
from .. import Utilities

class MetadataTF:
    """
Genotype metadata format
Metadata is a dataframe containing 'chromosome, position, id, allele_0, allele_1, allele_1_frequency'.
We are following https://mathgen.stats.ox.ac.uk/impute/1000GP%20Phase%203%20haplotypes%206%20October%202014.html,
we assume that 'allele_0' is the 'reference allele' and 'allele_1' is the alternate or effect allele.
The actual studies might have additional columns. i.e. In gtex id is not an rsid, but we have a n extra column with that if available.
    """
    CHROMOSOME=0
    POSITION=1
    ID=2
    ALLELE_0=3
    ALLELE_1=4
    ALLELE_1_FREQUENCY=5

    K_CHROMOSOME="chromosome"
    K_POSITION="position"
    K_ID="id"
    K_ALLELE_0="allele_0"
    K_ALLELE_1="allele_1"
    K_ALLELE_1_FREQUENCY="allele_1_frequency"

    order = [(CHROMOSOME, K_CHROMOSOME), (POSITION, K_POSITION), (ID, K_ID), (ALLELE_0, K_ALLELE_0), (ALLELE_1, K_ALLELE_1), (ALLELE_1_FREQUENCY, K_ALLELE_1_FREQUENCY)]


class MetadataTFE(MetadataTF):
    RSID=6
    K_RSID = "rsid"
    order = MetadataTF.order + [(RSID, K_RSID)]

class Genotype():
    """
Abstraction for holding genotype from a study.
Variant is a list, and each entry is a numpy array with the variant's individual values
Metadata is a dataframe in the format of 'MetadataTF' class; columns might be null/NA
    """
    def __init__(self, variants, metadata):
        self.variants = variants
        self.metadata = metadata

    def get_variants_metadata(self, variants=None):
        return _get_variants_metadata(self.metadata, variants)

    def get_variants(self, variants=None, to_pandas=True, specific_individuals=None):
        return _get_variants(self.metadata, self.variants, variants, to_pandas, specific_individuals)


def _metadata_from_raw(data, extra_columns=None):
    columns = [key for order,key in MetadataTF.order]
    if extra_columns: columns = columns + extra_columns
    return Utilities.to_dataframe(data, columns)

def _to_minor_allele_frequency(genotype):
    """
Use at your own risk
    """
    g_ = copy.deepcopy(genotype)
    m_ = g_.metadata
    clause_ = m_.allele_1_frequency > 0.5
    F = MetadataTF
    m_.loc[clause_,[F.K_ALLELE_0, F.K_ALLELE_1]] = m_.loc[clause_,[F.K_ALLELE_1, F.K_ALLELE_0]].values
    m_.loc[clause_, F.K_ALLELE_1_FREQUENCY] = 1 - m_.loc[clause_, F.K_ALLELE_1_FREQUENCY].values
    variants = g_.variants

    for i,swap in enumerate(clause_):
        if swap:
            variants[i] = 2 - variants[i]
    return g_

def _get_variants_metadata(m, variants=None):
    if variants is not None:
        m = m.loc[m.id.isin(variants)]
        m = m.assign(order = pandas.Categorical(m.id, categories=variants, ordered=True))
        m = m.sort_values(by="order")
        m = m.drop("order", axis=1)
    return m

def _to_variants(metadata, variants, to_pandas):
    if metadata.shape[0] != len(variants):
        raise RuntimeError("We expect metadata and variants of the same length (Actually the same order)")
    d = {x.id:variants[i] for i,x in enumerate(metadata.itertuples())}
    if to_pandas:
        d = pandas.DataFrame(d)
        d = d[metadata.id.values]
    return d

def _get_variants(metadata, variants, variant_list, to_pandas, specific_individuals):
    if not variant_list:
        return _to_variants(metadata, variants, to_pandas)
    if specific_individuals is not None:
        raise RuntimeError("Not implemented at the moment")

    index = {}
    _v = {x for x in variant_list}
    for i, tuple in enumerate(metadata.itertuples()):
        if tuple.id in _v: index[tuple.id] = i

    variant_list = [x for x in variant_list if x in index]
    v = [variants[index[x]] for x in variant_list]
    m = _get_variants_metadata(metadata, variant_list)
    return _to_variants(m, v, to_pandas)

def _monoallelic_by_frequency(metadata):
    #Ugly code to avoid pandas grouping inefficiency
    key = []
    for e in metadata.itertuples():
        key.append("{}_{}".format(e.chromosome, e.position))
    metadata["key"] = key

    v = {}
    _v = set()
    for e in metadata.itertuples():
        _i = e.key
        if _i in v:
            _v.add(_i)
        else:
            v[_i] = []
        v[_i].append((_i, e.id, e.allele_1_frequency if e.allele_1_frequency<0.5 else 1-e.allele_1_frequency))

    remove=set()
    for id in _v:
        entries = v[id]
        largest = entries[0]

        for i in range(1, len(entries)):
            entry = entries[i]
            if entry[2] > largest[2]:
                largest = entry

        for e in entries:
            if e[1] != largest[1]:
                remove.add(e[1])

    d = metadata[~metadata.id.isin(remove)]
    d.drop("key", axis=1, inplace=True)
    return d
