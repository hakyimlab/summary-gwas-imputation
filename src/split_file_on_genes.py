#!/usr/bin/env python
__author__ = "alvaro barbeira"
import os
import re
import logging
import gzip
import pandas
from genomic_tools_lib import Utilities, Logging
from genomic_tools_lib.miscellaneous import PandasHelpers

KEY="split"
PATH="path"

def get_split_map(annotation, splits, format):
    split_map = {}
    chromosomes = ["chr{}".format(i) for i in range(1,23)]
    for c in chromosomes:
        d = annotation.loc[annotation.chromosome == c]
        for split in range(0, splits):
            s = PandasHelpers.sub_batch(d, splits, split)
            for t in s.itertuples():
                split_map[t.gene_id] = {}
                key = "{}_{}".format(c, split)
                split_map[t.gene_id][KEY] = key
                split_map[t.gene_id][PATH] = format.format(c, split)
    return split_map

intron_re=re.compile("chr(\d+):(\d+):(\d+):(.*):(.*)")
def convert_to_intron(k):
    s = intron_re.search(k)
    return "intron_{}_{}_{}".format(s.group(1), s.group(2), s.group(3))

def run(args):
    logging.info("Loading annotation")
    annotation = pandas.read_table(args.gene_annotation, usecols=["chromosome", "gene_id"])

    logging.info("Creating split map")
    split_map = get_split_map(annotation, args.splits, args.output_format)

    conversion = None
    if args.key_conversion == "INTRON":
        conversion = convert_to_intron

    Utilities.ensure_requisite_folders(args.output_format)
    logging.info("Processing")
    last_key=None
    last_file=None
    wrote_header=set()
    unmmapped=set()
    with gzip.open(args.input) as f:
        header = f.readline()
        for i,line in enumerate(f):
            comps = line.decode().split()
            gene_id = comps[0]
            if conversion:
                gene_id = conversion(gene_id)

            if not gene_id in split_map:
                if not gene_id in unmmapped:
                    logging.log(9, "Unmapped gene %s", gene_id)
                    unmmapped.add(gene_id)
                continue

            _split = split_map[gene_id]
            if last_key != _split[KEY]:
                if last_file:
                    logging.log(9, "Closing %s", last_key)
                    last_file.close()
                last_key = _split[KEY]
                logging.log(9, "Opening %s", last_key)
                last_file = gzip.open(_split[PATH], "a")
                if not last_key in wrote_header:
                    last_file.write(header)
                    wrote_header.add(last_key)
            last_file.write(line)

    logging.info("Finished processing")

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser("Convert gencode file into a table file")
    parser.add_argument("-gene_annotation", help="Where to load file from")
    parser.add_argument("-output_format", help="Where to save")
    parser.add_argument("-input", help="input file to split")
    parser.add_argument("--key_conversion", help="optional, [INTRON, NONE]")
    parser.add_argument("-splits", help="Specify multiple key-value pairs to specify format conversion", type=int)
    parser.add_argument("-parsimony", help="Logging verbosity (actually loquacity)", type=int, default=10)
    args = parser.parse_args()

    Logging.configure_logging(args.parsimony)
    run(args)