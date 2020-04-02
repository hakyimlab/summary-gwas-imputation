#!/usr/bin/env python3
__author__ = "alvaro barbeira"
import logging
import os
import sqlite3
import pandas

from genomic_tools_lib import Logging, Utilities
from genomic_tools_lib.data_management import KeyedDataSource
from genomic_tools_lib.miscellaneous import Models
from genomic_tools_lib.file_formats import Miscellaneous

def run(args):
    if os.path.exists(args.output):
        logging.info("Output exists. Nope.")
        return

    Utilities.ensure_requisite_folders(args.output)

    logging.info("Loading variant annotation")
    variants = KeyedDataSource.load_data(args.variant_annotation, "variant_id", args.rsid_column)

    logging.info("Loading data annotation")
    if len(args.data_annotation) == 1:
        data_annotation = pandas.read_table(args.data_annotation[0])
        data_annotation = data_annotation[["gene_id", "gene_name", "feature_type", "gene_type"]][data_annotation.feature_type == "gene"].drop_duplicates()
    elif len(args.data_annotation) == 2:
        data_annotation = pandas.read_table(args.data_annotation[0])
        data_annotation = data_annotation[["gene_id", "gene_name", "feature_type", "gene_type"]][
        data_annotation.feature_type == args.data_annotation[1]].drop_duplicates()
    else:
        raise  RuntimeError("Unsupported annotation length")

    logging.info("Loading model_input")
    data = pandas.read_table(args.model_input, usecols=["gene_id", "gene_name", "variant", "weight"])

    logging.info("Processing")
    if args.model_filter and args.model_filter[1] == "PIP":
        w = Miscellaneous.dapg_signals(args.model_filter[0], float(args.model_filter[2]), variants)
        w = w.rename(columns={"gene":"gene_id", "variant_id":"variant"})
        data = data.merge(w[["gene_id", "variant"]], on=["gene_id", "variant"])

    v = pandas.DataFrame([(k,variants[k]) for k in data.variant.drop_duplicates()], columns=["variant", "rsid"])
    v.loc[v.rsid == ".", "rsid"] = v.loc[v.rsid == ".", "variant"]
    weights = data.merge(v, on="variant")
    weights = weights.assign(
        ref_allele = weights.variant.str.replace("(.*)_(.*)_(.*)_(.*)_b38", lambda x: x.group(3)),
        eff_allele=weights.variant.str.replace("(.*)_(.*)_(.*)_(.*)_b38", lambda x: x.group(4)))
    weights = weights.rename(columns={"variant":"varID", "gene_id":"gene"})[["gene", "rsid", "varID", "ref_allele", "eff_allele", "weight"]]

    extra = data.groupby("gene_id").size().to_frame("n.snps.in.model").reset_index()
    extra = extra.merge(data_annotation[["gene_id", "gene_name", "gene_type"]], on="gene_id")
    extra["pred.perf.pval"] = None
    extra["pred.perf.qval"] = None
    extra["pred.perf.R2"] = None
    extra = extra[["gene_id", "gene_name", "gene_type", "n.snps.in.model", "pred.perf.R2", "pred.perf.pval", "pred.perf.qval"]].rename(columns={"gene_id":"gene", "gene_name":"genename"})

    logging.info("Saving db")
    Models.create_model_db(args.output, extra, weights)

    logging.info("Done")

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser("Build predition model from data")
    parser.add_argument("-data_annotation", nargs="+")
    parser.add_argument("-variant_annotation")
    parser.add_argument("-model_input")
    parser.add_argument("--model_filter", nargs="+")
    parser.add_argument("-output")
    parser.add_argument("-parsimony", type=int, default=logging.INFO)
    parser.add_argument("-rsid_column", help = "Column name with rsid variant identifiers", default = "rs_id_dbSNP150_GRCh38p7")
    args = parser.parse_args()

    Logging.configure_logging(args.parsimony)
    run(args)
