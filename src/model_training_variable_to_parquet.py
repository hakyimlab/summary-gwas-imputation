#!/usr/bin/env python
__author__ = "alvaro barbeira"

import logging

from timeit import default_timer as timer

from genomic_tools_lib import Logging
from genomic_tools_lib import Utilities
from genomic_tools_lib.file_formats import ModelTraining, Parquet

def run(args):
    start = timer()
    Utilities.ensure_requisite_folders(args.parquet_output)
    logging.info("Loading variable")
    variables = ModelTraining.load_variable_file(args.variable_file)
    logging.info("Saving")
    Parquet.save_variable(args.parquet_output, variables)
    end = timer()
    logging.info("Finished in %s", str(end-start))

if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser("Convert model training format data to parquet format ")
    parser.add_argument("-variable_file", help="Folder where genotype files are")
    parser.add_argument("-parquet_output", help="Parquet file to save")
    parser.add_argument("-parsimony", help="Log parsimony level. 1 is everything being logged. 10 is only high level messages, above 10 will hardly log anything", default = "10")

    args = parser.parse_args()

    Logging.configure_logging(int(args.parsimony))

    run(args)