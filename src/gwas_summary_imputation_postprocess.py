#!/usr/bin/env python
__author__ = "alvaro barbeira"
import os
import re
import logging
from timeit import default_timer as timer

import pandas
import numpy
from scipy import stats

from genomic_tools_lib import Logging, Utilities
from genomic_tools_lib.miscellaneous import Genomics

COLUMN_ORDER = ["variant_id", "panel_variant_id", "chromosome", "position", "effect_allele", "non_effect_allele",
                # Original_build (TODO),
                "current_build", "frequency", "sample_size", "zscore", "pvalue", "effect_size", "standard_error",
                "imputation_status", "n_cases"]


def gwas_k(d):
    r = []
    for t in d.itertuples():
        _r = None
        try:
            _r = "{}_{}".format(t.chromosome, int(t.position))
        except Exception as e:
            logging.log(6, "Error for {}_{}".format(t.chromosome, t.position))
        r.append(_r)
    return r


def process_original_gwas(args, imputed):
    logging.info("Processing GWAS file %s", args.gwas_file)
    g = pandas.read_table(args.gwas_file)
    g = g.assign(current_build="hg38", imputation_status="original")[COLUMN_ORDER]
    # Remember the palindromic snps are to be excluded from the input GWAS;
    logging.info("Read %d variants", g.shape[0])

    if not args.keep_all_observed:
        if args.keep_criteria == "GTEX_VARIANT_ID":
            g = g.loc[~ g.panel_variant_id.isin(imputed.panel_variant_id)]
        elif args.keep_criteria == "CHR_POS":
            g = g.assign(k = gwas_k(g))
            imputed = imputed.assign(k = gwas_k(imputed))
            g = g.loc[~ g.k.isin({x for x in imputed.k})]
            g.drop("k", axis=1, inplace=True)
            imputed.drop("k", axis=1, inplace=True)
        else:
            raise RuntimeError("Unsupported keep option")
        logging.info("Kept %d variants as observed", g.shape[0])

    g = pandas.concat([g, imputed])[COLUMN_ORDER]
    logging.info("%d variants", g.shape[0])

    logging.info("Filling median")
    g = Genomics.fill_column_to_median(g, "sample_size", numpy.int32)

    logging.info("Sorting by chromosome-position")
    g = Genomics.sort(g)

    logging.info("Saving")
    Utilities.save_dataframe(g, args.output)

    return g[["panel_variant_id"]]


def process_imputed(args):
    r = re.compile(args.pattern)
    files = sorted([x for x in os.listdir(args.folder) if r.search(x)])
    result =[]
    for i, file in enumerate(files):
        logging.info("Processing imputed %s", file)
        p = os.path.join(args.folder, file)
        g = pandas.read_table(p)
        if g.shape[0] == 0:
            logging.info("Empty set of results for %s", p)
            continue

        g.drop(["n", "n_indep", "most_extreme_z"], axis=1, inplace=True)
        g.rename(columns={"effect_allele_frequency": "frequency", "status": "imputation_status"}, inplace=True)
        g = g.assign(pvalue=2 * stats.norm.sf(numpy.abs(g.zscore)), effect_size=numpy.nan, standard_error=numpy.nan,
                     sample_size=numpy.nan, current_build="hg38", n_cases=numpy.nan)
        g = g[COLUMN_ORDER]
        result.append(g)
    result = pandas.concat(result)
    logging.info("Processed %d imputed variants", result.shape[0])
    return result


def run(args):
    if os.path.exists(args.output):
        logging.info("Output exists. Nope.")
        return

    start = timer()
    logging.info("Beginning process")

    imputed = process_imputed(args)
    process_original_gwas(args, imputed)

    end = timer()
    logging.info("Finished in %s seconds", str(end - start))


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser("Post summary imputation results")
    parser.add_argument("-gwas_file",
                        help="GWAS file. For the moment, uniform-formatted hg38-based files are accepted.")
    parser.add_argument("-folder", help="How far to extend in each direction when searching for variants")
    parser.add_argument("-pattern", help="Work only with one chromosome")
    parser.add_argument("-output", help="Where to save stuff")
    parser.add_argument("--keep_all_observed", help="If an imputed gwas is present in the observed values",
                        action="store_true")
    parser.add_argument("--keep_criteria",
                        help="Discard original entries according to match by: CHR_POS or GTEX_VARIANT_ID",
                        default="GTEX_VARIANT_ID", type=str)
    parser.add_argument("-parsimony",
                        help="Log verbosity level. 1 is everything being logged. 10 is only high level messages, above 10 will hardly log anything",
                        default="10")

    args = parser.parse_args()

    Logging.configure_logging(int(args.parsimony))

    run(args)
