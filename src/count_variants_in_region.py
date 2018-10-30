#!/usr/bin/env python
import genomic_tools_lib.miscellaneous.Genomics

__author__ = "alvaro barbeira"
import logging

from pyarrow import parquet as pq

from genomic_tools_lib.miscellaneous import Genomics
from genomic_tools_lib import Logging, Utilities

def count_variants(chromosome, start, end, vf, m, last_chromosome, args):
    try:
        chromosome = int(chromosome.split("chr")[1])
        start = int(start)
        end = int(end)
        if chromosome != last_chromosome:
            logging.info("Reading chromosome %d", chromosome)
            m = vf.read_row_group(chromosome-1).to_pandas()
            last_chromosome = chromosome
            if args.frequency_filter:
                logging.log(9, "Filtering by frequency")
                m = m[(m.allele_1_frequency > args.frequency_filter) & (m.allele_1_frequency < 1-args.frequency_filter)]
        v = Genomics.entries_for_window(chromosome, start, end, m)
        count = v.shape[0]
    except:
        count = "NA"
    return count, m, last_chromosome

def run(args):
    logging.info("Starting process")

    vf = pq.ParquetFile(args.parquet_genotype_metadata)
    m = None
    last_chromosome = None

    r = []
    for i, line in Utilities.iterate_file(args.regions):
        if i==0: continue
        comps = line.strip().split()
        count, m, last_chromosome = count_variants(comps[0], comps[1], comps[2], vf, m, last_chromosome, args)
        r.append((comps[0], comps[1], comps[2], count))

    r = Utilities.to_dataframe(r, ["chromosome", "start", "end", "count"])
    Utilities.save_dataframe(r, args.output)
    logging.info("Finished process")

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser("Generate BSLMM runs on study")
    parser.add_argument("-regions", help="Ld rediongs.")
    parser.add_argument("-parquet_genotype_metadata", help="Parquet Genotype variant metadata file")
    parser.add_argument("-window", help="How far to extend in each direction when searching for variants", type=int)
    parser.add_argument("-frequency_filter", help="Skip variants with frequency (below f) or above (1-f)", type=float)
    parser.add_argument("-output")
    parser.add_argument("-parsimony", help="Log verbosity level. 1 is everything being logged. 10 is only high level messages, above 10 will hardly log anything", default = "10")

    args = parser.parse_args()

    Logging.configure_logging(int(args.parsimony))

    run(args)

