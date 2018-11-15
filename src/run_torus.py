#!/usr/bin/env python
__author__ = "alvaro barbeira"

import os
import logging
from timeit import default_timer as timer

from genomic_tools_lib import Logging
from genomic_tools_lib.external_tools.torus import RunTorus, Utilities as RunTorusUtilities

def run(args):
    start = timer()
    if os.path.exists(args.intermediate_folder):
        logging.info("Intermediate folder exists. Nope.")
        return

    if os.path.exists(args.output_folder):
        logging.info("Output folder exists. Nope.")
        return

    context = RunTorusUtilities.context_from_args(args)
    RunTorus.run_torus(context)
    end = timer()
    logging.info("Ran torus in %s", str(end-start))

if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser("Convert model training format data to parquet format ")
    parser.add_argument("-torus_command", help="Which torus to run")
    parser.add_argument("-frequency_filter", help="If provided, will ignore those entries that don't satisfy a frequency criteria", type=float)
    parser.add_argument("-snp_annotation_file", help="file with snp annotation")
    parser.add_argument("-gene_annotation", help="file with gene annotation. Add optional string -parsed- or -gencode- (default) to specify format", nargs="+")
    parser.add_argument("-eqtl", help="File with eqtl")
    parser.add_argument("-eqtl_mode", help="eQTL processing mode: -SQTL- or -EQTL-(default)", default="EQTL")
    parser.add_argument("-intermediate_folder", help="Folder where scratch data will be written to")
    parser.add_argument("-output_folder", help="Folder where torus output will be written to. Must not be a subfolder of intermediate.")
    parser.add_argument("-keep_intermediate_folder", help="Wether to keep torus intermediate stuff", action="store_true")
    parser.add_argument("-snp_annotation_from_parquet_metadata", help="Load a genotype study metadata (parquet format) for the required information")
    parser.add_argument("-parsimony", help="Log verbosity level. 1 is everything being logged. 10 is only high level messages, above 10 will hardly log anything", default = "10")

    args = parser.parse_args()

    Logging.configure_logging(int(args.parsimony))

    run(args)