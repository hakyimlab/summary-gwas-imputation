#!/usr/bin/env python
__author__ = "alvaro barbeira"
import logging
import pandas
import gzip
from genomic_tools_lib import Utilities, Logging

def get_gene_whitelist(args):
    g = pandas.read_table(args.gencode)
    return {x for x in g.gene_id.values}

def run(args):
    logging.info("Acquiring whitelist")
    #whitelist = get_gene_whitelist(args)

    logging.info("Processing...")
    Utilities.ensure_requisite_folders(args.output)
    with gzip.open(args.output, "w") as _o:
        _o.write("gene\tcluster\tchromosome\tstart\tend\n".encode())
        with gzip.open(args.input) as _i:
            for i,line in enumerate(_i):
                if i==0: continue
                comps = line.decode().strip().split()
                d = comps[3].split(":")
                o = "{}\t{}\t{}\t{}\t{}\n".format(d[4], d[3], d[0], d[1], d[2]).encode()
                _o.write(o)

    logging.info("Finished")

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser("Patch gene name column")
    parser.add_argument("-input", help="Where to load file from")
    parser.add_argument("-output", help="Where to save")
    parser.add_argument("-gencode", help="Gencode file to filter which genes to keep")
    parser.add_argument("-parsimony", help="Logging parsimony", type=int, default=10)
    args = parser.parse_args()

    Logging.configure_logging(args.parsimony)
    run(args)