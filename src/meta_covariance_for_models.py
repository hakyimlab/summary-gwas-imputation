__author__ = "alvaro barbeira"

import logging
import os
import re
import sqlite3
import pandas
import numpy
import gzip

from pyarrow import parquet as pq

from genomic_tools_lib import Logging, Utilities
from genomic_tools_lib.miscellaneous import matrices, PandasHelpers
from genomic_tools_lib.file_formats import Parquet

def get_file_map(args):
    r = re.compile(args.parquet_genotype_pattern)
    files = os.listdir(args.parquet_genotype_folder)
    files = {int(r.search(f).groups()[0]):os.path.join(args.parquet_genotype_folder, f) for f in files if r.search(f)}
    p = {}
    for k,v in files.items():
        logging.log(9, "opening:%s", v)
        g = pq.ParquetFile(v)
        p[k] = g
    return p

n_ = re.compile("^(\d+)$")

def get_gene_variant_list(model_db_folder, pattern):
    logic = Utilities.file_logic(model_db_folder, pattern)
    g = []
    for i,l in enumerate(logic.itertuples()):
        logging.log(9, "Opening db %d/%d: %s ", i+1, logic.shape[0], l.name)
        with sqlite3.connect(l.path) as connection:
            w = pandas.read_sql("select * from weights;", connection)[["gene", "varID", "rsid"]]
            g.append(w)
    g = pandas.concat(g).drop_duplicates()
    return g

def run(args):
    if os.path.exists(args.output):
        logging.info("Output already exists, either delete it or move it")
        return

    logging.info("Getting parquet genotypes")
    file_map = get_file_map(args)

    logging.info("Getting variants")
    gene_variants = get_gene_variant_list(args.model_db_folder, args.model_db_file_pattern)
    genes = list(gene_variants.gene.drop_duplicates())

    Utilities.ensure_requisite_folders(args.output)

    logging.info("Processing")
    with gzip.open(args.output, "w") as f:
        f.write("GENE RSID1 RSID2 VALUE\n".encode())
        for i,g in enumerate(gene_variants.gene.drop_duplicates()):
            logging.log(9, "Proccessing %i/%i:%s", i+1, len(genes), g)
            w = gene_variants[gene_variants.gene == g]
            chr_ = w.varID.values[0].split("_")[0].split("chr")[1]
            if not n_.search(chr_):
                logging.log(9, "Unsupported chromosome: %s", chr_)
                continue

            dosage = file_map[int(chr_)]
            d = Parquet._read(dosage, columns=w.varID.values, skip_individuals=True)
            var_ids = list(d.keys())
            if args.output_rsids:
                ids = [x for x in pandas.DataFrame({"varID": var_ids}).merge(w[["varID", "rsid"]], on="varID").rsid.values]
            else:
                ids = var_ids
            c = numpy.cov([d[x] for x in var_ids])
            c = matrices._flatten_matrix_data([(w.gene.values[0], ids, c)])
            for entry in c:
                l = "{} {} {} {}\n".format(entry[0], entry[1], entry[2], entry[3])
                f.write(l.encode())
    logging.info("Finished building covariance.")

if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser("Build S-MultiXcan covariance from models")
    parser.add_argument("-parquet_genotype_folder", help="Parquet Genotype folder")
    parser.add_argument("-parquet_genotype_pattern", help="Pattern to detect parquet genotypes by chromosome")
    parser.add_argument("-model_db_folder", help="Folder containing models")
    parser.add_argument("-model_db_file_pattern", help="Regexp to parse file names")
    parser.add_argument("-output", help="Where to save stuff")
    parser.add_argument("--output_rsids", action="store_true")
    parser.add_argument("-parsimony", help="Log verbosity level. 1 is everything being logged. 10 is only high level messages, above 10 will hardly log anything", default = "10")

    args = parser.parse_args()

    Logging.configure_logging(int(args.parsimony))

    run(args)