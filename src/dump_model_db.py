import gzip
import pandas
import numpy
import logging
import re
import os

from genomic_tools_lib import Utilities, Logging
from genomic_tools_lib.miscellaneous import Models

__author__ = "alvaro barbeira"

def _read(folder, pattern):
    files = sorted([x for x in os.listdir(folder) if pattern.search(x)])
    files = [os.path.join(folder, x) for x in files]
    r =[]
    for file in files:
        try:
            r.append(pandas.read_table(file))
        except:
            logging.info("issue opening file %s", file)
    return pandas.concat(r)

def _read_2(input_prefix, stem):
    path_ = os.path.split(input_prefix)
    r = re.compile(path_[1] + stem)
    return _read(path_[0], r)

def run(args):

    logging.info("Loading models")
    weights, extra = Models.read_model(args.input)

    Utilities.save_dataframe(weights, args.output_prefix + "_weights.txt.gz")
    Utilities.save_dataframe(extra, args.output_prefix + "_extra.txt.gz")
    logging.info("Done")

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser("Dump")
    parser.add_argument("-input")
    parser.add_argument("-output_prefix")
    parser.add_argument("-parsimony", type=int, default=logging.INFO)

    args = parser.parse_args()
    Logging.configure_logging(args.parsimony)
    run(args)
