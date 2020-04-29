__author__ = "alvaro barbeira"
import logging
import pandas
import gzip
import re
from . import PandasHelpers

def discard_gtex_palindromic_variants(d):
    return d[~(d.id.str.contains("_A_T_") | d.id.str.contains("_T_A_") | d.id.str.contains("_C_G_") | d.id.str.contains("_G_C_"))]

def allele_key(d):
    def _a(x):
        if x.effect_allele < x.non_effect_allele:
            a1, a2 = x.effect_allele, x.non_effect_allele
        else:
            a1, a2 = x.non_effect_allele, x.effect_allele
        return "{}_{}".format(a1, a2)

    return d.apply(_a, axis=1)

def _genomic_index(d):
    i={}
    for k,e in enumerate(d.itertuples()):
        if not e.chromosome in i: i[e.chromosome] = {}
        _i = i[e.chromosome]

        if not e.position in _i:
            _i[e.position] = {}

        _i[e.position][e.non_effect_allele] = {}
        _i[e.position][e.non_effect_allele][e.effect_allele] = e.panel_variant_id

    return i

def _build_alignment(source, i):
    flip = []
    id = []
    for e in source.itertuples():
        _f = False
        _id = None
        if e.chromosome in i:
            _i = i[e.chromosome]
            if e.position in _i:
                _i = _i[e.position]
                if  e.effect_allele in _i and e.non_effect_allele in _i[e.effect_allele]:
                    _f = True
                    _id = _i[e.effect_allele][e.non_effect_allele]
                elif e.non_effect_allele in _i and e.effect_allele in _i[e.non_effect_allele]:
                    _id = _i[e.non_effect_allele][e.effect_allele]
                else:
                    pass
        flip.append(_f)
        id.append(_id)
    return flip, id

def to_number(series):
    s=[]
    for x in series:
        s_ = None
        try:
            s_ = float(x)
        except:
            pass
        s.append(s_)
    return pandas.Series(s, dtype=float)

def to_int(series):
    s=[]
    for x in series:
        s_ = "NA"
        try:
            s_ = int(x)
        except:
            pass
        s.append(s_)
    return pandas.Series(s)

def match(source, reference):
    #Ugly code to go faster than pandas's engine
    logging.log(9, "Indexing reference")
    i = _genomic_index(reference)
    logging.log(9, "building flip")
    flip, id = _build_alignment(source, i)

    aligned = source.assign(panel_variant_id = id)
    aligned.loc[flip, ["non_effect_allele", "effect_allele"]] = aligned.loc[flip, ["effect_allele", "non_effect_allele"]].values
    aligned.loc[flip, "zscore"] = - aligned.loc[flip, "zscore"]

    if "effect_size" in aligned.columns.values:
        aligned.loc[flip, "effect_size"] = - aligned.loc[flip, "effect_size"]

    if "frequency" in aligned.columns.values:
        aligned.loc[flip, "frequency"] = 1.0 - aligned.loc[flip, "frequency"]

    return aligned

def sort(d):
    chr_re_ = re.compile("chr(\d+)$")
    chr = [int(x.split("chr")[1]) if chr_re_.search(x) else None for x in d.chromosome]
    d = d.assign(chr = chr)
    d = d.sort_values(by=["chr", "position"])
    d.drop("chr", axis=1, inplace=True)
    return d

def fill_column_to_median(d, column, dtype=None):
    c = d[column]
    if c[c.isna()].shape[0] == c.shape[0]:
        logging.info("Warning: column %s has no defined values", column)
    else:
        m = c[~c.isna()].median()
        d.loc[c.isna(), column] = m
        if dtype:
            d[column] = d[column].astype(dtype)
    return d

def entries_for_window(chromosome, window_start, window_end, metadata):
    if window_start < 0:
        window_start = 0
    m = metadata[(metadata.position >= window_start) & (metadata.position <= window_end) & (metadata.chromosome == chromosome)]
    return m

def entries_for_gene_annotation(annotation, window, metadata):
    min_ = annotation.start - window
    max_ = annotation.end + window
    m = entries_for_window(annotation.chromosome, min_, max_, metadata)
    return m

def entries_for_split(chromosome, splits, split, metadata):
    metadata = metadata[metadata.chromosome == chromosome]
    return PandasHelpers.sub_batch(metadata, splits, split)
