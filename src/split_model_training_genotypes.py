#!/usr/bin/env python
__author__ = "alvaro barbeira"
import logging
import os
import gzip
import numpy
import pandas

from genomic_tools_lib import Utilities, Logging
from genomic_tools_lib.data_management import TextFileTools
from genomic_tools_lib.miscellaneous import Genomics

def build_regions(annotation, chromosome, sub_jobs, window):
    results=[]
    genes=[]
    for i in range(0, sub_jobs):
        s = Genomics.entries_for_split(chromosome, sub_jobs, i, annotation)
        start = numpy.min(s.start)- window
        if start<0: start=0
        end = numpy.max(s.end) + window
        results.append((i+1, start, end))
        genes.append(s.gene_id.values)
    return pandas.DataFrame(results, columns=["split", "start", "end"]), genes

def run(args):
    logging.info("Loading annotation")
    annotation = pandas.read_table(args.input_annotation)

    logging.info("Loading region")
    regions,genes = build_regions(annotation, args.chromosome, args.sub_jobs, args.window)

    file_name = os.path.split(args.input_file)[1]
    name = file_name.split(".txt.gz")[0]

    logging.info("Saving gene lists")
    gene_outputs = [os.path.join(args.output_folder, name) + "_{}_genes.txt.gz".format(i) for i in range(1, args.sub_jobs + 1)]
    for i,p in enumerate(gene_outputs):
        with gzip.open(p, "w") as f:
            genes_ = genes[i]
            for gene in genes_:
                f.write("{}\n".format(gene).encode())

    logging.info("Processing file")
    outputs = [os.path.join(args.output_folder,name)+"_{}.txt.gz".format(i) for i in range(1, args.sub_jobs+1)]

    Utilities.ensure_requisite_folders(outputs[0])
    output_files = [gzip.open(x, "w") for x in outputs]

    with gzip.open(args.input_file) as input_file:
        header = input_file.readline()
        for f in output_files:
            f.write(header)

        for i,line in enumerate(input_file):
            comps = line.decode().strip().split()
            pos = int(comps[0].split("_")[1])
            targets = regions[(regions.start<=pos)&(pos<regions.end)]
            for target in targets.itertuples():
                f = output_files[target.Index]
                f.write(line)

    logging.info("Finalizing output files")
    for f in output_files:
        f.close()

    logging.info("Finished")

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser("Convert GTEx splicing phenotype file into a table file")
    parser.add_argument("-input_annotation", help="Where to load annotation")
    parser.add_argument("-input_file", help="Where to load file from")
    parser.add_argument("-output_folder", help="Where to save")
    parser.add_argument("-chromosome", help="Where to load annotation", type=int)
    parser.add_argument("-sub_jobs", help="How much to split, basically", type=int)
    parser.add_argument("-window", help="Where to load file from", type=int, default=1000000)
    parser.add_argument("-verbosity", help="Logging verbosity (actually loquacity)", type=int, default=10)
    args = parser.parse_args()

    Logging.configure_logging(args.verbosity)
    run(args)