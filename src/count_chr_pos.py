__author__ = "alvaro barbeira"
import os
import logging
import pandas
import gzip
import csv

from genomic_tools_lib import Utilities, Logging


def duplicated_entries_in_file(path):
    r = {}
    with gzip.open(path) as f:
        f.readline()
        for i,l in enumerate(f):
            comps = l.decode().strip().split()
            key = "{}\t{}".format(comps[2], comps[3])
            if not key in r:
                r[key] = 0
            r[key] += 1
    r = {k:v for k,v in r.items() if v>1}
    return pandas.DataFrame({"key":list(r.keys()), "count": list(r.values())})

def duplicated_entries(folder):
    r = []
    files = sorted([x for x in os.listdir(folder) if "txt.gz" in x])
    for f in files:
        logging.log(9, "Processing %s", f)
        path = os.path.join(folder, f)
        r_ = duplicated_entries_in_file(path)
        r_["file"] = f
        if r_ is not None:
            r.append(r_)
            if r_.shape[0] > 0:
                logging.log(9, "Duplicated found!")
    return pandas.concat(r)

def run(args):
    d = duplicated_entries(args.input_folder)
    Utilities.ensure_requisite_folders(args.output)
    Utilities.save_dataframe(d, args.output, sep=",", quoting=csv.QUOTE_NONE)

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser("Check imputed gwas in a folder for duplicated chr-pos entries")
    parser.add_argument("-input_folder")
    parser.add_argument("-output")
    parser.add_argument("-parsimony", default=10, type=int)
    args = parser.parse_args()
    Logging.configure_logging(args.parsimony)

    run(args)