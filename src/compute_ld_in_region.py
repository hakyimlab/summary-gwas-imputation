#!/usr/bin/env python
__author__ = "alvaro barbeira"

import logging
import os
import re
import pandas
import numpy
import gzip
from timeit import default_timer as timer

from pyarrow import parquet as pq

from genomic_tools_lib import Logging, Utilities
from genomic_tools_lib.data_management import TextFileTools
from genomic_tools_lib.miscellaneous import matrices, PandasHelpers
from genomic_tools_lib.miscellaneous import Genomics
from genomic_tools_lib.file_formats import Parquet

class Context:
    def __init__(self, args):
        self.args = args
        self.file_map = None
        self.vmf = None
        self.of = None
        self.regions = None

    def get_genotype_file(self, chromosome):
        logging.info("Opening genotype for chromosome %d", chromosome)
        g = pq.ParquetFile(self.file_map[chromosome])
        return g

    def __enter__(self):
        logging.info("initializing resources")

        logging.info("Loading regions")
        regions = load_regions(self.args.region_file, self.args.chromosome)
        if args.sub_batches and args.sub_batch is not None:
            logging.log(9, "Selecting target regions from sub-batches")
            regions = PandasHelpers.sub_batch(regions, args.sub_batches, args.sub_batch)
        self.regions = regions

        logging.info("Opening variants metadata")
        self.vmf = pq.ParquetFile(args.parquet_genotype_metadata)

        logging.info("Creating destination")
        if args.text_output:
            if os.path.exists(args.text_output):
                raise  RuntimeError("Output exists. Nope.")
            Utilities.ensure_requisite_folders(args.text_output)
            self.of = TextFileTools.TextDataSink(args.text_output, [("region", "id1", "id2", "value")])
            self.of.initialize()
        elif args.text_output_folder:
            Utilities.maybe_create_folder(args.text_output_folder)
        else:
            raise RuntimeError("Unrecognized output specification")

        if (args.parquet_genotype_folder and args.parquet_genotype_pattern):
            self.file_map = get_file_map(args)
        else:
            raise RuntimeError("Unrecognized genotype specification")

        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        logging.info("finalizing resources")
        if self.of:
            self.of.finalize()

    def sink(self, cov, ids, region):
        logging.log(9, "Serializing covariance")
        _region = "{}_{}_{}_{}".format(region.name, region.chr, region.start, region.stop)
        if args.text_output:
            if args.dapg_output:
                raise RuntimeError("Not supported for this option")
            else:
                cov = matrices._flatten_matrix_data([(_region, ids, cov)])
                self.of.sink(cov)
        elif args.text_output_folder:
            if args.dapg_output:
                f = os.path.join(args.text_output_folder, _region) + ".txt.gz"
                with gzip.open(f, "w") as o:
                    for i in range(0, cov.shape[0]):
                        l = "\t".join(["{:0.4f}".format(x) for x in cov[i]]) + "\n"
                        o.write(l.encode())
                id = os.path.join(args.text_output_folder, _region) + ".id.txt.gz"
                with gzip.open(id, "w") as o:
                    l = "\n".join(ids).encode()
                    o.write(l)

            else:
                cov = matrices._flatten_matrix_data_2(ids, cov)
                cov = pandas.DataFrame(cov)[["id1", "id2", "value"]]
                f = os.path.join(args.text_output_folder, _region) + ".txt.gz"
                Utilities.save_dataframe(cov, f)

def get_file_map(args):
    r = re.compile(args.parquet_genotype_pattern)
    files = os.listdir(args.parquet_genotype_folder)
    files = {int(r.search(f).groups()[0]):os.path.join(args.parquet_genotype_folder, f) for f in files if r.search(f)}
    return files

def filter_by_frequency(vm, frequency):
    return vm.loc[(frequency < vm.allele_1_frequency) &
                  (vm.allele_1_frequency < 1 - frequency)]

def load_regions(path, chromosome):
    regions = pandas.read_table(path)
    regions = regions.assign(name = ["region_{}".format(x) for x in regions.index.values])
    regions.dropna(inplace=True)
    regions = regions.assign(start = regions.start.astype(numpy.int32), stop = regions.stop.astype(numpy.int32))
    if chromosome:
        regions = regions.loc[regions.chr == "chr{}".format(chromosome)]
    return regions

#TODO remove debug code
def by_chromosome(context, chromosome):
    vm = context.vmf.read_row_group(chromosome - 1).to_pandas()
    if args.frequency_filter:
        vm = filter_by_frequency(vm, args.frequency_filter)

    g = context.get_genotype_file(chromosome)

    regions = context.regions
    regions = regions[regions.chr == "chr{}".format(chromosome)]

    for i,region in enumerate(regions.itertuples()):
        logging.log(9, "Processing region in chr %d: %d/%d", chromosome, i+1, regions.shape[0])
        vmw = Genomics.entries_for_window(chromosome, region.start - args.window, region.stop + args.window, vm)
        ids = vmw.id.values
        logging.log(9, "%d variants", len(ids))
        d = Parquet._read(g, columns=ids, skip_individuals=True)
        d = numpy.array([d[x] for x in ids], dtype=numpy.float32)
        if context.args.standardise_geno:
            cov = numpy.corrcoef(d, ddof=1).astype(numpy.float32, copy=False)
        else:
            cov = numpy.cov(d).astype(numpy.float32, copy=False)
        logging.log(9, "%d rows", cov.shape[0])
        context.sink(cov, ids, region)


def run(args):
    start = timer()
    logging.info("Starting")

    with Context(args) as context:
        if args.chromosome:
            by_chromosome(context, args.chromosome)
        else:
            for chromosome in range(1,23):
                by_chromosome(context, chromosome)

    end = timer()
    logging.info("Ended in %s", str(end-start))

if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser("Generate BSLMM runs on study")
    parser.add_argument("-region_file", help="LD-independent regions.")
    parser.add_argument("-parquet_genotype", help="Parquet Genotype folder")
    parser.add_argument("-parquet_genotype_folder", help="Parquet Genotype folder")
    parser.add_argument("-parquet_genotype_pattern", help="Pattern to detect parquet genotypes by chromosome")
    parser.add_argument("-parquet_genotype_metadata", help="Parquet Genotype variant metadata file")
    parser.add_argument("-window", help="How far to extend in each direction when searching for variants", type=int, default=0)
    parser.add_argument("-chromosome", help="Work only with one chromosome", type=int)
    parser.add_argument("-text_output", help="Where to save stuff")
    parser.add_argument("-text_output_folder", help="Where to save stuff")
    parser.add_argument("--frequency_filter", help="Skip variants with frequency (below f) or above (1-f)", type=float)
    parser.add_argument("-sub_batches", help="Split the data into subsets", type=int)
    parser.add_argument("-sub_batch", help="only do this subset", type=int)
    parser.add_argument("-parsimony", help="Log verbosity level. 1 is everything being logged. 10 is only high level messages, above 10 will hardly log anything", default = "10")
    parser.add_argument("--dapg_output", help="Output matrices in DAP-G format", action="store_true")
    parser.add_argument("--standardise_geno", help="Standardise geno, or get correlation matrix", action="store_true")

    args = parser.parse_args()

    Logging.configure_logging(int(args.parsimony))

    run(args)