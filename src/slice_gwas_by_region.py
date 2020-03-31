#!/usr/bin/env python
__author__ = "alvaro barbeira"
import os
import logging

import numpy
import pandas

from genomic_tools_lib import Utilities, Logging

def run(args):
    if os.path.exists(args.output):
        logging.info("%s exists. Nope.", args.output)
        return

    Utilities.ensure_requisite_folders(args.output)

    logging.info("Loading regions")
    regions = pandas.read_table(args.region_file).rename(columns={"chr":"chromosome"})
    regions.dropna(inplace=True)
    regions.start = regions.start.astype(int)
    regions.stop = regions.stop.astype(int)

    logging.info("Loading gwas")
    gwas = pandas.read_table(args.gwas_file, usecols=["panel_variant_id", "chromosome", "position", "zscore"])
    gwas.dropna(inplace=True)

    logging.info("Processing")
    sliced = []
    for i,region in enumerate(regions.itertuples()):
        logging.log(8, "Processing region %d", i+1)
        if numpy.isnan(region.start) or numpy.isnan(region.stop) or \
                (type(region.chromosome) != str and numpy.isnan(region.chromosome)):
            logging.log(8, "skipping incomplete region")
            continue
        slice = gwas[(gwas.chromosome == region.chromosome) & (gwas.position >= region.start) & (gwas.position < region.stop)]
        slice = slice.sort_values(by = "position")
        if slice.shape[0] == 0:
            continue
        slice = slice.assign(region = "region-{}-{}-{}".format(region.chromosome, region.start, region.stop), r=i)

        slice = slice[["panel_variant_id", "region", "r", "zscore"]]
        sliced.append(slice)

    sliced = pandas.concat(sliced).sort_values(by="r")
    if args.output_format == "dapg":
        sliced.region = sliced.r.apply(lambda x: "region{}".format(x))
        sliced = sliced.drop(["r"], axis=1)
        Utilities.save_dataframe(sliced, args.output, header=False)
    elif args.output_format == "gtex_eqtl":
        sliced = sliced.assign(gene_id = sliced.region, variant_id=sliced.panel_variant_id, tss_distance = numpy.nan, ma_samples = numpy.nan, ma_count= numpy.nan, maf = numpy.nan, pval_nominal = numpy.nan, slope= sliced.zscore, slope_se=1)
        sliced = sliced[["gene_id", "variant_id", "tss_distance", "ma_samples", "ma_count", "maf", "pval_nominal", "slope", "slope_se"]]
        Utilities.save_dataframe(sliced, args.output, header=True)
    logging.info("Finished slicing gwas")

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser("Group gwas results by region")
    parser.add_argument("-region_file", help="Non-overlapping regions")
    parser.add_argument("-gwas_file", help="GWAS file, in fixed format (imputed) for now")
    parser.add_argument("-output", help="Where to save the result")
    parser.add_argument("-parsimony", help="How much logging to output", type=int, default=10)
    parser.add_argument("--output_format", default="dapg")

    args = parser.parse_args()
    Logging.configure_logging(args.parsimony)
    run(args)