#!/usr/bin/env python
__author__ = "alvaro barbeira"
import os
import re
import logging
from timeit import default_timer as timer

import pandas
import numpy
from scipy import stats

from genomic_tools_lib import  Logging, Utilities
from genomic_tools_lib.file_formats.gwas import Utilities as GWASUtilities

COLUMN_ORDER=["variant_id", "panel_variant_id", "chromosome", "position", "effect_allele", "non_effect_allele", #Original_build (TODO),
               "current_build", "frequency", "sample_size", "zscore", "pvalue", "effect_size", "standard_error", "imputation_status"]

def _gwas_k(x):
    r = None
    try:
        r = "{}_{}".format(x.chromosome, int(x.position))
    except Exception as e:
        logging.log(6, "Error for {}_{}".format(x.chromosome, x.position))
    return r

def process_original_gwas(args, imputed_keys):
    logging.info("Processing GWAS file %s", args.gwas_file)
    g = pandas.read_table(args.gwas_file)
    #Remember the palindromic snps are to be excluded from the input GWAS;
    logging.info("Read %d variants", g.shape[0])
    if not args.keep_all_observed:
        if args.keep_criteria == "GTEX_VARIANT_ID":
            g = g.loc[~ g.panel_variant_id.isin(imputed_keys)]
        elif args.keep_criteria == "CHR_POS":
            g["k"] = g.apply(_gwas_k, axis=1)
            g = g.loc[~ g.k.isin(imputed_keys)]
            g.drop("k", axis=1, inplace=True)
        else:
            raise RuntimeError("Unsupported keep option")
        logging.info("Kept %d variants as observed", g.shape[0])
    g = g.assign(current_build = "hg38", imputation_status="original")[COLUMN_ORDER]
    Utilities.save_dataframe(g, args.output, mode="a", header=False)

    return g[["panel_variant_id"]]

def process_imputed(args):
    r = re.compile(args.pattern)
    files = sorted([x for x in os.listdir(args.folder) if r.search(x)])
    count = 0
    keys = set()
    for i,file in enumerate(files):
        logging.info("Processing imputed %s", file)
        p = os.path.join(args.folder, file)
        g = pandas.read_table(p)
        if g.shape[0] == 0:
            logging.info("Empty set of results for %s", p)
            continue
        count += g.shape[0]

        #Fast dropping of observed values
        #g = g.merge(observed_ids, on="panel_variant_id", how="left", copy=False, indicator=True)
        #g = g[g._merge == "left_only"]

        g.drop(["n", "n_indep", "most_extreme_z"], axis=1, inplace=True)
        g.rename(columns={"effect_allele_frequency": "frequency", "status":"imputation_status"}, inplace=True)
        g = g.assign(pvalue = 2*stats.norm.sf(numpy.abs(g.zscore)), effect_size= numpy.nan, standard_error = numpy.nan, sample_size=numpy.nan, current_build = "hg38")
        g = g[COLUMN_ORDER]
        Utilities.save_dataframe(g, args.output, mode="a" if i>0 else "w", header = i==0)
        if not args.keep_all_observed:
            if args.keep_criteria == "GTEX_VARIANT_ID":
                keys.update(g.panel_variant_id.values)
            elif args.keep_criteria == "CHR_POS":
                chr_pos = g.apply(lambda x: "{}_{}".format(x.chromosome, int(x.position)),axis=1)
                keys.update(chr_pos)
            else:
                raise RuntimeError("Unsupported keep option")


    logging.info("Processed %d imputed variants", count)
    return keys

def run(args):
    if os.path.exists(args.output):
        logging.info("Output exists. Nope.")
        return

    start = timer()
    logging.info("Beginning process")

    imputed_keys = process_imputed(args)
    process_original_gwas(args, imputed_keys)

    end = timer()
    logging.info("Finished in %s seconds", str(end-start))


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser("Collect summary imputation results")
    parser.add_argument("-gwas_file", help="GWAS file. For the moment, uniform-formatted hg38-based files are accepted.")
    parser.add_argument("-folder", help="How far to extend in each direction when searching for variants")
    parser.add_argument("-pattern", help="Work only with one chromosome")
    parser.add_argument("-output", help="Where to save stuff")
    parser.add_argument("--keep_all_observed", help="If an imputed gwas is present in the observed values", action="store_true")
    parser.add_argument("--keep_criteria", help="Discard original entries according to match by: CHR_POS or GTEX_VARIANT_ID", default="CHR_POS", type=str)
    parser.add_argument("-parsimony", help="Log verbosity level. 1 is everything being logged. 10 is only high level messages, above 10 will hardly log anything", default = "10")

    args = parser.parse_args()

    Logging.configure_logging(int(args.parsimony))

    run(args)