#!/usr/bin/env python
__author__ = "alvaro barbeira"

import os
import logging

from timeit import default_timer as timer

from genomic_tools_lib import Logging
from genomic_tools_lib.external_tools.dap import Utilities as DAPUtilities, RunDAP
from genomic_tools_lib import Utilities

def run(args):
    start = timer()

    if os.path.exists(args.output_folder):
        logging.info("Output folder exists. Nope.")
        return

    if os.path.exists(args.intermediate_folder):
        logging.info("Intermediate folder exists. Nope.")
        return

    stats = []

    context = DAPUtilities.context_from_args(args)
    available_genes = context.get_available_genes()

    for i,gene in enumerate(available_genes):
        if args.MAX_M and i==args.MAX_M:
            break
        _start = timer()
        logging.log(8, "Processing %i/%i:%s", i+1, len(available_genes), gene)
        _stats = RunDAP.run_dap(context, gene)
        _end = timer()
        logging.log(7, "Elapsed: %s", str(_end - _start))
        stats.append(_stats)

    end = timer()
    logging.info("Ran DAP in %s seconds" % (str(end - start)))

    Utilities.ensure_requisite_folders(args.output_folder)
    stats_ = args.stats_name if args.stats_name else "stats.txt"
    stats_path = os.path.join(args.output_folder, stats_)
    stats = RunDAP.data_frame_from_stats(stats).fillna("NA")
    Utilities.save_dataframe(stats, stats_path)

if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser("Generate DAP runs on study")
    parser.add_argument("-dap_command", help="Which gemma command to run")
    parser.add_argument("-frequency_filter", help="If provided, restrict to variants satisfying MAF criteria", type=float)
    parser.add_argument("-grid_file", help="Which grid file to use")
    parser.add_argument("-priors_folder", help="Folder containing (torus) priors")
    parser.add_argument("-options", help="DAP-G command line options", action="append", nargs=2)
    parser.add_argument("-intermediate_folder", help="Folder to use as scratch space")
    parser.add_argument("-gene_annotation", help="File describing genes")
    parser.add_argument("-parquet_genotype", help="Parquet Genotype file")
    parser.add_argument("-parquet_genotype_metadata", help="Parquet Genotype variant metadata file")
    parser.add_argument("-parquet_phenotype", help="Parquet phenotypes")
    parser.add_argument("-parquet_covariate", help= "Parquet covariates")
    parser.add_argument("-window", help="How far to extend in each direction when searching for variants", type=int)
    parser.add_argument("-output_folder", help="Where will the model output weights be saved")
    parser.add_argument("-chromosome", help="Split the data into subsets", type=int)
    parser.add_argument("-sub_batches", help="Split the data into subsets", type=int)
    parser.add_argument("-sub_batch", help="only do this subset", type=int)
    parser.add_argument("--keep_intermediate_folder", help="don't delete the intermediate stuff", action='store_true')
    parser.add_argument("--MAX_M", help="Run for up o this many genes", type=int, default=None)
    parser.add_argument("--stats_name", help="name for run stats")
    parser.add_argument("-parsimony", help="Log verbosity level. 1 is everything being logged. 10 is only high level messages, above 10 will hardly log anything", default = "10")

    args = parser.parse_args()

    Logging.configure_logging(int(args.parsimony))

    run(args)