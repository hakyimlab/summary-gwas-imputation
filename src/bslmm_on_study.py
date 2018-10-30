#!/usr/bin/env python
__author__ = "alvaro barbeira"

import os
import shutil
import logging
import pandas

from timeit import default_timer as timer

from genomic_tools_lib import Logging, Utilities
from genomic_tools_lib.external_tools.gemma import Utilities as GEMMAUtilities, RunGEMMA

def run(args):
    start = timer()
    if os.path.exists(args.intermediate_folder):
        logging.info("Intermediate folder already exists. Nope.")
        return
    else:
        os.makedirs(args.intermediate_folder)

    if not (args.output_stats or args.output_weights or args.output_covariance or args.output_hyperparameters):
        logging.info("Specify at least one output")
        return

    if args.output_weights: Utilities.ensure_requisite_folders(args.output_weights)
    if args.output_stats: Utilities.ensure_requisite_folders(args.output_stats)
    if args.output_covariance: Utilities.ensure_requisite_folders(args.output_covariance)
    if args.output_hyperparameters: Utilities.ensure_requisite_folders(args.output_hyperparameters)

    weights = []
    covariance = []
    stats = []
    hyperparameters = []
    context = GEMMAUtilities.context_from_args(args)
    n = len(context.get_available_genes())
    for i,gene in enumerate(context.get_available_genes()):
        start_ = timer()
        logging.log(8, "Processing %d/%d:%s", i+1,n, gene)
        weights_,  covariance_, hyperparameters_, stats_ = RunGEMMA.run_gemma(context, gene)

        end_ = timer()
        logging.log(7, "Elapsed: %s", str(end_-start_))
        if args.output_weights: weights.append(weights_)
        if args.output_covariance: covariance.append(covariance_)
        if args.output_hyperparameters: hyperparameters.append(hyperparameters_)
        if args.output_stats: stats.append(stats_)

    if args.output_weights:
        weights = pandas.concat(weights)
        Utilities.save_dataframe(weights, args.output_weights)

    if args.output_stats:
        stats = RunGEMMA.dataframe_from_stats(stats).fillna("NA")
        Utilities.save_dataframe(stats, args.output_stats)

    if args.output_covariance:
        covariance = RunGEMMA.dataframe_from_covariance_data(covariance).fillna("NA")
        Utilities.save_dataframe(covariance, args.output_covariance)

    if args.output_hyperparameters:
        hyperparameters = RunGEMMA.dataframe_from_hyperparameters(hyperparameters).fillna("NA")
        Utilities.save_dataframe(hyperparameters, args.output_hyperparameters)

    shutil.rmtree(args.intermediate_folder)

    end = timer()
    logging.info("Ran BSLMM in %s seconds" % (str(end - start)))

if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser("Generate BSLMM runs on study")
    parser.add_argument("-gene_annotation", help="File describing genes")
    parser.add_argument("-gemma_command", help="Which gemma command to run")
    parser.add_argument("-intermediate_folder", help="What gemma to run")
    parser.add_argument("-parquet_genotype", help="Parquet Genotype file")
    parser.add_argument("-parquet_genotype_metadata", help="Parquet Genotype variant metadata file")
    parser.add_argument("-parquet_phenotype", help="Parquet phenotypes")
    parser.add_argument("-parquet_covariate", help="Parquet covariates")
    parser.add_argument("-window", help="How far to extend in each direction when searching for variants", type=int)
    parser.add_argument("-output_weights", help="Where will the model output weights be saved")
    parser.add_argument("-output_covariance", help="Where will the model covariance be saved")
    parser.add_argument("-output_stats", help="Where will the model output stats be saved")
    parser.add_argument("-output_hyperparameters", help="Where will the model output hyperparameters be saved")
    parser.add_argument("-sub_batches", help="Split the data into subsets", type=int)
    parser.add_argument("-sub_batch", help="only do this subset", type=int)
    parser.add_argument("-verbosity", help="Log verbosity level. 1 is everything being logged. 10 is only high level messages, above 10 will hardly log anything", default = "10")

    args = parser.parse_args()

    Logging.configure_logging(int(args.verbosity))

    run(args)