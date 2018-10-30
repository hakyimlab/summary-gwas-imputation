#!/usr/bin/env python
__author__ = "alvaro barbeira"
import os
import logging

import pandas

from genomic_tools_lib import Utilities, Logging

def run(args):
    if os.path.exists(args.output):
        logging.info("%s exists. Nope.", args.output)
        return

    logging.info("Loading regions")
    regions = pandas.read_table(args.region_file).rename(columns={"chr":"chromosome"})
    regions.dropna(inplace=True)

    logging.info("Loading gwas")
    gwas = pandas.read_table(args.gwas_file, usecols=["panel_variant_id", "chromosome", "position", "zscore"])
    gwas.dropna(inplace=True)

    logging.info("Processing")
    sliced = []
    for i,region in enumerate(regions.itertuples()):
        logging.log(8, "Processing region %d", i+1)
        slice = gwas[(gwas.chromosome == region.chromosome) & (gwas.position >= region.start) & (gwas.position < region.stop)]
        slice = slice.sort_values(by = "position")
        if slice.shape[0] == 0:
            continue
        slice = slice.assign(region = "region{}".format(i+1), r=i)[["panel_variant_id", "region", "r", "zscore"]]
        sliced.append(slice)

    sliced = pandas.concat(sliced).sort_values(by="r").drop(["r"], axis=1)
    Utilities.save_dataframe(sliced, args.output, header=False)
    logging.info("Finished slicing gwas")

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser("Group gwas results by region")
    parser.add_argument("-region_file", help="Non-overlapping regions")
    parser.add_argument("-gwas_file", help="GWAS file, in fixed format (imputed) for now")
    parser.add_argument("-output", help="Where to save the result")
    parser.add_argument("-parsimony", help="How much logging to output", type=int, default=10)

    args = parser.parse_args()
    Logging.configure_logging(args.parsimony)
    run(args)