import gzip
import pandas
import sqlite3
import logging
import re
import os

from genomic_tools_lib import Utilities, Logging

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
    logging.info("Loading model summaries")
    extra = _read_2(args.input_prefix, "_(.*)_summary.txt.gz")
    extra = extra[(extra["pred.perf.pval"] < 0.05) & (extra.rho_avg >0.1)]
    logging.info("Saving extra")
    db = args.output_prefix + ".db"
    with sqlite3.connect(db) as conn:
        extra.to_sql("extra", conn, index=False)

    logging.info("Loading weights")
    weights = _read_2(args.input_prefix, "_(.*)_weights.txt.gz")
    weights = weights[weights.gene.isin(extra.gene)]
    logging.info("Saving weights")
    with sqlite3.connect(db) as conn:
        weights.to_sql("weights", conn, index=False)

    logging.info("Processing covariances")
    genes = {x for x in extra.gene}

    path_ = os.path.split(args.input_prefix)
    r = re.compile(path_[1] + "_(.*)_covariance.txt.gz")
    files = sorted([x for x in os.listdir(path_[0]) if r.search(x)])
    files = [os.path.join(path_[0], x) for x in files]
    cov = args.output_prefix + ".txt.gz"
    with gzip.open(cov, "w") as cov_:
        cov_.write("GENE RSID1 RSID2 VALUE\n".encode())
        for f in files:
            with gzip.open(f) as f_:
                f_.readline()
                for l in f_:
                    gene = l.decode().strip().split()[0]
                    if not gene in genes:
                        continue
                    cov_.write(l)
    logging.info("Done")


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser("Compile model database and covariance from model training output")
    parser.add_argument("-input_prefix")
    parser.add_argument("-output_prefix")
    parser.add_argument("-parsimony", type=int, default=logging.INFO)

    args = parser.parse_args()
    Logging.configure_logging(args.parsimony)
    run(args)