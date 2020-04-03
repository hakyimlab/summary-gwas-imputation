#!/usr/bin/env python
__author__ = "alvaro barbeira"

import os
import logging
from timeit import default_timer as timer

import numpy
import pandas
from scipy import stats

from genomic_tools_lib import Logging, Utilities
from genomic_tools_lib.summary_imputation import Utilities as SummaryImputationUtilities, SummaryInputation

def postprocess_results(results):
    results = results.rename(columns={"id":"variant_id", "status":"imputation_status"})
    results = results.assign(pvalue = 2*stats.norm.sf(numpy.abs(results.zscore)))
    # results = results.assign(standard_error = None).assign(effect_size=None).assign(sample_size=None).assign(current_build="hg38")
    # results = results[["variant_id", "panel_variant_id", "chromosome", "position",
    #                    "effect_allele", "non_effect_allele", "frequency", "pvalue", "zscore", "effect_size", "standard_error", "sample_size", "imputation_status"]]
    return results

def run_by_variant(args):
    context = SummaryImputationUtilities.context_from_args(args)

    r, s = [], []
    variant_metadata = context.get_target_variants_metadata()
    for i, variant in enumerate(variant_metadata.itertuples()):
        logging.log(8, "Variant %d/%d:%s", i+1, variant_metadata.shape[0], variant.id)
        _r, _s = SummaryInputation.gaussian(context, variant)
        r.append(_r)
        s.append(_s)

    results = SummaryInputation.dataframe_from_results(r,s)
    results = postprocess_results(results)

    return results

def run_by_region(args):
    context = SummaryImputationUtilities.context_by_region_from_args(args)

    results = []
    regions = context.get_target_regions()
    for i, region in enumerate(regions.itertuples()):
        logging.log(9, "Processing region %d/%d [{}, {}]".format(region.start, region.end), i+1, regions.shape[0])
        _r = SummaryInputation.gaussian_by_region(context, region)
        results.append(_r)

    results = pandas.concat(results)
    return results

def run(args):
    if os.path.exists(args.output):
        logging.info("Output exists. Nope.")
        return

    start = timer()
    logging.info("Beginning process")
    if args.by_region_file:
        results = run_by_region(args)
    else:
        results = run_by_variant(args)

    Utilities.ensure_requisite_folders(args.output)
    Utilities.save_dataframe(results, args.output)

    end = timer()
    logging.info("Finished in %s seconds", str(end-start))

if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser("Impute summary statistics.")
    parser.add_argument("-by_region_file", help="If provided, run imputation grouping by regions. Much faster.")
    parser.add_argument("-parquet_genotype", help="Parquet Genotype file")
    parser.add_argument("-parquet_genotype_metadata", help="Parquet Genotype variant metadata file")
    parser.add_argument("-gwas_file", help="GWAS file. For the moment, uniform-formatted hg38-based files are accepted.")
    parser.add_argument("-window", help="How far to extend in each direction when searching for variants", type=int, default=0)
    parser.add_argument("-chromosome", help="Work only with one chromosome", type=int)
    parser.add_argument("-output", help="Where to save stuff")
    parser.add_argument("-cutoff", help="naive cutoff when performing SVD", type=float, default=0)
    parser.add_argument("-regularization", help="Ridge-like regularization for matrix inversion", type=float, default=0)
    parser.add_argument("-frequency_filter", help="Skip variants with frequency (below f) or above (1-f)", type=float)
    parser.add_argument("-sub_batches", help="Split the data into subsets", type=int)
    parser.add_argument("-sub_batch", help="only do this subset", type=int)
    parser.add_argument("-containing", help="only do this subset", type=int)
    parser.add_argument("--keep_palindromic_imputation", help="Report imputed values istead of (original+flipped) values for palindromic snps", action="store_true")
    parser.add_argument("--use_palindromic_snps", help="Use palindromic variants when imputing", action="store_true")
    parser.add_argument("--standardise_dosages", help="Standardise dosages before computing (i.e. use correlation matrix)", action="store_true")
    parser.add_argument("--cache_variants", help="Save variants in memory instead of loading every time, when running by variant", action="store_true")
    parser.add_argument("-parsimony", help="Log verbosity level. 1 is everything being logged. 10 is only high level messages, above 10 will hardly log anything", default = "10")

    args = parser.parse_args()

    Logging.configure_logging(int(args.parsimony))

    run(args)
