import os
import logging
from timeit import default_timer as timer
from pyarrow import parquet as pq
import pandas
import numpy
import re
import sys
import rpy2.robjects as robjects
from rpy2.robjects import pandas2ri
import subprocess
from genomic_tools_lib import Logging
from genomic_tools_lib.file_formats import Parquet


def run(args):
    pass



if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser()

    parser.add_argument('-out_dir', help="directory for results", required=True)
    parser.add_argument('-region', help="filepath for regions", required=False)
    parser.add_argument('-chr', help="Chromosome number", type=int,
                        required=True)
    parser.add_argument('-annotation', help="Parquet annotation file",
                        required=True)
    parser.add_argument('-geno', help="Parquet genotype file")
    parser.add_argument('-pheno', help="Parquet phenotype file")
    parser.add_argument('--out_split_by', help="What should output be split by?"
                                               " Options are 'region' or "
                                               "'pheno'")
    parser.add_argument('--n_batches', type=int)
    parser.add_argument('--batch', type=int)
    parser.add_argument('--parsimony', type=int, default=7)
    parser.add_argument('--compress', help="Gzip all resulting files after writing",
                        default=False, action='store_true')
    parser.add_argument('--MAX_R', type=int)

    args = parser.parse_args()

    Logging.configure_logging(args.parsimony, target=sys.stdout, with_date=True)
    run(args)