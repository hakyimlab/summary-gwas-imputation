__author__ = "alvaro barbeira"
from collections import namedtuple

import logging
import math
import numpy
import pandas
from numpy import dot as d

from ..miscellaneous import Math
from .. import Utilities
from ..file_formats.gwas import Utilities as GWASUtilities
from ..individual_data.Utilities import _StudyBasedContext
from ..individual_data import Utilities as StudyUtilities
from ..miscellaneous import Genomics

class _Context(_StudyBasedContext):
    def get_gwas_slice(self, variants_metadata, variant): raise RuntimeError("Not implemented")
    def get_cutoff(self): raise RuntimeError("Not implemented")
    def get_regularization(self): raise RuntimeError("Not implemented")
    def get_freq_filter(self): raise RuntimeError("Not implemented")
    def standardise_dosages(self): raise RuntimeError("Not implemented")
    def keep_palindromic_imputation(self): raise RuntimeError("Not Implemented")
    def use_palindromic_snps(self): raise RuntimeError("Not Implemented")

class _VariantContext(_Context):
    def get_target_variants_metadata(self): raise RuntimeError("Not implemented")

class _RegionContext(_Context):
    def get_target_regions(self): raise RuntimeError("Not implemented")


Results = namedtuple("Results", ["variant_id", "panel_variant_id", "chromosome", "position", "effect_allele", "non_effect_allele", "frequency", "zscore", "variance", "imputation_status"])

AdditionalStats = namedtuple("AdditionalStats", ["n","n_indep","most_extreme_z"])

def _additional_stats(n="NA", n_indep="NA", most_extreme_z="NA"):
    return (n, n_indep, most_extreme_z)

def _result(variant, z, status):
    _r = (variant.rsid, variant.id, variant.chromosome, variant.position, variant.effect_allele, variant.non_effect_allele, variant.effect_allele_frequency, z, None, status)
    return _r

def _freq_filter(d, freq_filter):
    f = numpy.mean(d) / 2
    return f < freq_filter or (1 - freq_filter) < f

def _freq_filter_variant(variant, freq_filter):
    return variant.effect_allele_frequency < freq_filter or (1 - freq_filter) < variant.effect_allele_frequency

def _filter_variants_by_freq(v, f):
    v = v[(v.effect_allele_frequency > f) & (v.effect_allele_frequency < (1-f))]
    return v

#optimist function
def _get(variants, ids, cutoff, regularization, f=True):
    geno = [variants[x] for x in ids]
    cov = numpy.cov(geno)
    sigma = cov[:len(ids)-1,:len(ids)-1]
    rho = cov[-1:,:-1][0]
    sigma_inv, n_indep, eigen = Math.crpinv(sigma, cutoff, regularization)
    w = d(rho, sigma_inv)
    s = math.sqrt(d(rho,w)) if f else math.sqrt(d(w,d(sigma,w)))
    return w, s, n_indep, eigen

def _get_z(variants, ids, gwas_slice, cutoff, regularization, f=True):
    w, s, n_indep, eigen = _get(variants, ids, cutoff, regularization, f=f)
    z = d(w, gwas_slice.zscore) / s if s > 0 else None
    return  z, s, n_indep, eigen

def _get_variants(context, ids):
    variants = context.get_variants(ids)
    del variants["individual"]
    if context.standardise_dosages():
        #ordinary loop to save memory
        for k in variants.keys():
            variants[k] = Math.standardize(variants[k])
    return variants

########################################################################################################################
def _gaussian(context, variant):
    freq_filter = context.get_freq_filter()
    if freq_filter and _freq_filter_variant(variant, freq_filter):
        _r = _result(variant, numpy.nan, "freq_rejected")
        return _r, _additional_stats()

    variants_metadata = context.get_variants_metadata()
    window = context.get_window()

    variants_metadata_window = Genomics.entries_for_window(variant.chromosome, variant.position - window, variant.position + window, variants_metadata)
    if freq_filter:
        variants_metadata_window = _filter_variants_by_freq(variants_metadata_window, freq_filter)

    gwas_slice = context.get_gwas_slice(variants_metadata_window)
    if gwas_slice.shape[0] == 0:
        _r = _result(variant, numpy.nan, "no_gwas_gtex_intersection")
        return _r, _additional_stats()

    _existing = gwas_slice[(gwas_slice.chromosome == variant.chromosome) & (gwas_slice.position == variant.position)]
    if _existing.shape[0] > 0:
        logging.log(8, "Variant %s already exists", _existing.variant_id.values[0])
        gwas_slice = gwas_slice[gwas_slice.variant_id != _existing.variant_id.values[0]]

    _idx = numpy.argmax(numpy.abs(gwas_slice.zscore.values))
    most_extreme_z = gwas_slice.zscore[_idx]

    ids = list(gwas_slice.panel_variant_id.values) + [variant.id]

    variants = _get_variants(context, ids)

    if len({x for x in variants[variant.id]}) == 1:
        _r = _result(variant, numpy.nan, "bad_variant")
        return _r, _additional_stats()

    #Finally do the imputation
    z, s, n_indep, eigen = _get_z(variants, ids, gwas_slice, context.get_cutoff(), context.get_regularization())
    if s == 0:
        _r = _result(variant, numpy.nan, "bad_variance")
        return _r, _additional_stats(len(eigen), n_indep, most_extreme_z)

    _r = _result(variant, z, "imputed")

    return _r, _additional_stats(len(eigen), n_indep, most_extreme_z)

def gaussian(context, variant):
    try:
        _r, _a = _gaussian(context, variant)
    except Exception as e:
        _r, _a = _result(variant, None, "unexpected_error"), _additional_stats()
    return _r, _a

########################################################################################################################

def _error_region(context, region):
    variants_metadata = context.get_variants_metadata()
    variants_metadata_window = Genomics.entries_for_window(region.chromosome, region.start, region.end, variants_metadata)
    r = []
    s = []
    for v in variants_metadata_window.itertuples():
        r.append(_result(v, None, "bad_region"))
        s.append(_additional_stats())

    return dataframe_from_results(r,s)

# def _results(use_variants_metadata_window, variants_metadata_window, zscores, most_extreme_z, n_indep, error_status):
#     results = []
#     stats = []
#
#     if use_variants_metadata_window is not None and zscores is not None:
#         for i,variant in enumerate(use_variants_metadata_window.itertuples()):
#             _r = _result(variant, zscores[i], "imputed")
#             results.append(_r)
#             _s = _additional_stats(use_variants_metadata_window.shape[0], n_indep, most_extreme_z)
#             stats.append(_s)
#
#     error_status = error_status if error_status else "bad_variant"
#     rest = variants_metadata_window[~variants_metadata_window.id.isin({x for x in use_variants_metadata_window.id})] if use_variants_metadata_window is not None else variants_metadata_window
#     for variant in rest.itertuples():
#         _r = _result(variant, numpy.nan, error_status)
#         results.append(_r)
#         _s = _additional_stats()
#         stats.append(_s)
#
#     return results, stats


def _get_multi(geno, typed, cutoff, regularization):
    cov = numpy.cov(geno)
    sigma_tt = cov[:typed.shape[0],:typed.shape[0]]
    sigma_it = cov[typed.shape[0]:,:typed.shape[0]]
    sigma_inv, n_indep, eigen = Math.crpinv(sigma_tt, cutoff, regularization)
    w = d(sigma_it, sigma_inv)
    zscore = d(w, typed.zscore)

    _w = d(sigma_it, sigma_inv)
    variance = numpy.sum(numpy.multiply(sigma_it,_w), axis=1)
    return zscore, variance, sigma_tt.shape[0], n_indep

def _trim_to_region(d, region):
    return d[(d.position >= region.start) & (d.position < region.end)]

def _post_pocess(context, untyped, full_gwas_slice, gwas_slice, region):
    if context.keep_palindromic_imputation():
        pass
    elif context.use_palindromic_snps():
        pass
    else:
        logging.log(8, "Palindromic Post processing")
        #TODO
        #In the following, maybe it is the same to allign to gtex alleles. Check.
        full_gwas_slice = full_gwas_slice.assign(chromosome = "chr"+full_gwas_slice.chromosome.astype(str))
        keys =["chromosome", "position", "effect_allele", "non_effect_allele"]

        #So... we assume that we'll always have a panel id, other way it's unusable
        measured_z = full_gwas_slice[keys + ["panel_variant_id", "zscore"]]
        measured_z = measured_z.assign(abs_z = numpy.abs(measured_z.zscore)).drop(columns=["zscore"])
        unambiguous_measured_z = GWASUtilities.discard_ambiguous(measured_z)
        ambiguous_measured_z = measured_z.loc[~measured_z.panel_variant_id.isin(unambiguous_measured_z.panel_variant_id)]

        untyped = untyped.assign(sign = numpy.sign(untyped.zscore))
        key = untyped[keys +["sign"]]
        key = key.merge(ambiguous_measured_z, on=keys)
        key = key.assign(measured_z = key.abs_z * key.sign)[["panel_variant_id", "measured_z"]]
        untyped = untyped.merge(key, on="panel_variant_id", how="left")
        untyped.loc[~untyped.measured_z.isnull(), "zscore"] = untyped.loc[~untyped.measured_z.isnull(), "measured_z"]
        #untyped.loc[~untyped.measured_z.isnull(), "imputation_status"] = "original"
        untyped = untyped.drop(columns=["measured_z", "sign"])

    untyped = untyped[list(Results._fields+AdditionalStats._fields)]
    return untyped


def _gaussian_by_region(context, region):
    logging.log(8, "Roll out imputation")
    freq_filter = context.get_freq_filter()
    variants_metadata = context.get_variants_metadata()
    window = context.get_window()

    variants_metadata_window = Genomics.entries_for_window(region.chromosome, region.start - window, region.end + window, variants_metadata)
    use_variants_metadata_window = _filter_variants_by_freq(variants_metadata_window, freq_filter) if freq_filter else variants_metadata_window

    # top = use_variants_metadata_window.groupby(["chromosome", "position"]).\
    #     apply(lambda x: x.sort_values(by="effect_allele_frequency", ascending=False)).\
    #     reset_index(drop=True).groupby(["chromosome", "position"]).head(1)
    # use_variants_metadata_window = use_variants_metadata_window[use_variants_metadata_window.id.isin(top.id)]
    # debug
    # k = use_variants_metadata_window[["chromosome", "position", "id"]].groupby(["chromosome", "position"]).aggregate(["count"])
    # k = k.reset_index()
    # k.columns = k.columns.droplevel(1)


    full_gwas_slice = context.get_gwas_slice(use_variants_metadata_window)
    if context.use_palindromic_snps():
        gwas_slice = full_gwas_slice
    else:
        gwas_slice = GWASUtilities.discard_ambiguous(full_gwas_slice)

    if gwas_slice.shape[0] == 0:
        logging.info("No GWAS intersection for region")
        return dataframe_from_results([], [])

    _idx = numpy.argmax(numpy.abs(gwas_slice.zscore.values))
    most_extreme_z = gwas_slice.zscore.values[_idx]

    logging.log(8, "Preparing data")
    typed = use_variants_metadata_window[use_variants_metadata_window.id.isin({x for x in gwas_slice.panel_variant_id})]
    typed = pandas.merge(typed, gwas_slice[["panel_variant_id", "zscore"]].rename(columns={"panel_variant_id": "id"}))
    untyped = use_variants_metadata_window[~use_variants_metadata_window.id.isin({x for x in gwas_slice.panel_variant_id})]
    use_variants_metadata_window = pandas.concat([typed, untyped], sort=False)

    ids = [x for x in use_variants_metadata_window.id.values]
    variants = _get_variants(context, ids)
    geno = [variants[x] for x in ids]

    logging.log(8, "Performing imputation")
    zscore, variance, n, n_indep = _get_multi(geno, typed, context.get_cutoff(), context.get_regularization())

    untyped = untyped.rename(columns={"id":"panel_variant_id"}).\
        assign(zscore=zscore, variance=variance, imputation_status="imputed", n=n, n_indep=n_indep, most_extreme_z=most_extreme_z, chromosome="chr"+untyped.chromosome.astype(str)).\
        rename(columns={"rsid":"variant_id", "effect_allele_frequency":"frequency"})[list(Results._fields+AdditionalStats._fields)]

    # The snps on the border of the region go to the start of the next
    untyped = _trim_to_region(untyped, region)
    logging.log(8, "Postprocessing")
    imputed = _post_pocess(context, untyped, full_gwas_slice, gwas_slice, region)

    return imputed

def gaussian_by_region(context, region):
    results = dataframe_from_results([], [])
    try:
        results = _gaussian_by_region(context, region)
    except Exception as e:
        logging.info("Error for region ({},{},{}): {}".format(region.chromosome, region.start, region.end, repr(e)))
    return  results


########################################################################################################################
def dataframe_from_results(r, s):
    d = Utilities.to_dataframe(r, list(Results._fields))
    a = Utilities.to_dataframe(s, list(AdditionalStats._fields))
    return pandas.concat([d,a], axis=1)
