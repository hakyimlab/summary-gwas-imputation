#!/usr/bin/env python
__author__ = "alvaro barbeira"
import pandas
import os
import logging
from timeit import default_timer as timer

from genomic_tools_lib import Logging, Utilities
from genomic_tools_lib.external_tools.coloc import Coloc

def run(args):
    Coloc.initialize(args.coloc_script)
    if os.path.exists(args.output):
        logging.info("Output exists. Nope.")
        return
    start = timer()

    logging.info("Loading gwas")
    gwas = Coloc.read_gwas(args.gwas, args.gwas_sample_size, args.gwas_mode)

    streamer = Coloc.eqtl_streamer(args.eqtl, gwas)

    results = []
    logging.info("Beggining process")
    MAX_N=args.MAX_N
    for i, d in enumerate(streamer):
        gene = d.gene_id.values[0]
        logging.log(9, "Processing gene %s", gene)
        eqtl = Coloc.get_eqtl(d, args.eqtl_sample_size, args.eqtl_mode)
        r = Coloc.coloc_on_gwas_eqtl(gene, gwas, eqtl, args.gwas_mode, args.eqtl_mode, args.p1, args.p2, args.p12)
        results.append(r)
        if MAX_N and i > MAX_N:
            logging.info("Early exit")
            break

    logging.info("Saving")
    results = Coloc.results_to_dataframe(results)
    Utilities.ensure_requisite_folders(args.output)
    Utilities.save_dataframe(results, args.output)
    end = timer()
    logging.info("Finished COLOC in %s seconds" % (str(end - start)))

if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser("Run COLOC")
    parser.add_argument("-coloc_script", help="Optional override of R's coloc infrastructure")
    parser.add_argument("-gwas", help="Which gwas to run")
    parser.add_argument("-gwas_mode", help="options in [pvalue, bse, zscore_1]", default="bse")
    parser.add_argument("-eqtl", help="Which eQTL to run")
    parser.add_argument("-eqtl_mode", help="options in [pvalue, bse, zscore_1]", default="bse")
    parser.add_argument("-gwas_sample_size", help="either 'FROM_GWAS' (default) or integer sample size", default="FROM_GWAS")
    parser.add_argument("-eqtl_sample_size", help="eQTL number of samples", type=int)
    parser.add_argument("-p1", type=float, default=1e-4)
    parser.add_argument("-p2", type=float, default=1e-4)
    parser.add_argument("-p12", type=float, default=1e-5)
    parser.add_argument("-output", help="Folder where torus output will be written to. Must not be a subfolder of intermediate.")
    parser.add_argument("-keep_intermediate_folder", help="Wether to keep torus intermediate stuff", action="store_true")
    parser.add_argument("-snp_annotation_from_parquet_metadata", help="Load a genotype study metadata (parquet format) for the required information")
    parser.add_argument("-parsimony", help="Log verbosity level. 1 is everything being logged. 10 is only high level messages, above 10 will hardly log anything", default = "10")
    parser.add_argument("-MAX_N", type=int)

    args = parser.parse_args()

    Logging.configure_logging(int(args.parsimony))

    run(args)