#!/usr/bin/env python
__author__ = "alvaro barbeira"

import logging
import pandas
from genomic_tools_lib import Utilities, Logging
from genomic_tools_lib.data_management import KeyedDataSource

def run(args):
    r_ = pandas.read_csv if ".csv" in args.input else pandas.read_table
    sep = "," if ".csv" in args.output else "\t"

    logging.info("Loading gene table")
    g = KeyedDataSource.load_data(args.gene_table, "gene_id", "gene_name")

    logging.info("Loading input")
    i = r_(args.input)

    gene_name = []
    for t in i.itertuples():
        gene_name.append(g[t.gene])
    i["gene_name"] = gene_name

    logging.info("saving")
    Utilities.save_dataframe(i, args.output, sep=sep)

    logging.info("Done")

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser("Patch gene name column")
    parser.add_argument("-input", help="Where to load file from")
    parser.add_argument("-output", help="Where to save")
    parser.add_argument("-gene_table", help="Which types of genes to keep")
    parser.add_argument("-verbosity", help="Logging verbosity (actually loquacity)", type=int, default=10)
    args = parser.parse_args()

    Logging.configure_logging(args.verbosity)
    run(args)