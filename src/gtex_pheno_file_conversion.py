#!/usr/bin/env python
__author__ = "alvaro barbeira"
import logging
import os
import gzip

from genomic_tools_lib import Utilities, Logging
from genomic_tools_lib.data_management import TextFileTools

def _to_line(comps, index, row_name):
    return "\t".join([row_name] + [x for i, x in enumerate(comps) if i in index])

def input_generator(input_file, samples):
    introns = set()
    for i,line in Utilities.iterate_file(input_file):
        comps = line.strip().split()
        if i==0:
            _index = [i for i,x in enumerate(comps) if x in samples]
            yield _to_line(comps, _index, "NAME")
            continue

        name = comps[3]
        name_ = name.split(":")
        name = "_".join(["intron", name_[0].split("chr")[1], name_[1], name_[2]])
        gene = name_[4]
        if name in introns:
            continue
        introns.add(name)

        yield _to_line(comps, _index, name)


def run(args):
    if os.path.exists(args.output):
        logging.info("Output exists. Nope.")
        return

    logging.info("Loading samples")
    samples = {x for x in TextFileTools.load_list(args.samples_whitelist)}

    logging.info("Processing file")
    Utilities.ensure_requisite_folders(args.output)
    Utilities.write_iterable_to_file(input_generator(args.input_file, samples), args.output)

    logging.info("Finished")



if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser("Convert GTEx splicing phenotype file into a table file")
    parser.add_argument("-input_file", help="Where to load file from")
    parser.add_argument("-samples_whitelist", help="Which types of genes to keep")
    parser.add_argument("-output", help="Where to save")
    parser.add_argument("-verbosity", help="Logging verbosity (actually loquacity)", type=int, default=10)
    args = parser.parse_args()

    Logging.configure_logging(args.verbosity)
    run(args)