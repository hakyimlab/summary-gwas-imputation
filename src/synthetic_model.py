#!/usr/bin/env python3
__author__ = "alvaro barbeira"
import logging
import os
import sqlite3
import pandas

from genomic_tools_lib import Logging, Utilities
from genomic_tools_lib.data_management import KeyedDataSource
from genomic_tools_lib.miscellaneous import Models

def run(args):
    if os.path.exists(args.output):
        logging.info("Output exists. Nope.")
        return

    Utilities.ensure_requisite_folders(args.output)

    logging.info("Loading variant annotation")
    variants = KeyedDataSource.load_data(args.variant_annotation, "variant_id", "rs_id_dbSNP150_GRCh38p7")

    logging.info("Loading data annotation")
    data_annotation = pandas.read_table(args.data_annotation)
    data_annotation = data_annotation[["gene_id", "gene_name", "feature_type", "gene_type"]][data_annotation.feature_type == "gene"].drop_duplicates()

    logging.info("Loading model_input")
    data = pandas.read_table(args.model_input, usecols=["gene_id", "gene_name", "variant", "weight"])

    extra = data.groupby("gene_id").size().to_frame("n.snps.in.model").reset_index()
    extra = extra.merge(data_annotation[["gene_id", "gene_name", "gene_type"]], on="gene_id")
    extra["pred.perf.pval"] = None
    extra["pred.perf.qval"] = None
    extra["pred.perf.R2"] = None
    extra = extra[["gene_id", "gene_name", "gene_type", "n.snps.in.model", "pred.perf.R2", "pred.perf.pval", "pred.perf.qval"]].rename(columns={"gene_id":"gene", "gene_name":"genename"})

    v = pandas.DataFrame([(k,variants[k]) for k in data.variant.drop_duplicates()], columns=["variant", "rsid"])
    v.loc[v.rsid == ".", "rsid"] = v.loc[v.rsid == ".", "variant"]
    weights = data.merge(v, on="variant")
    weights = weights.assign(
        ref_allele = weights.variant.str.replace("(.*)_(.*)_(.*)_(.*)_b38", lambda x: x.group(3)),
        eff_allele=weights.variant.str.replace("(.*)_(.*)_(.*)_(.*)_b38", lambda x: x.group(4)))
    weights = weights.rename(columns={"variant":"varID", "gene_id":"gene"})[["gene", "rsid", "varID", "ref_allele", "eff_allele", "weight"]]

    logging.info("Saving db")
    Models.create_model_db(args.output, extra, weights)

    logging.info("Done")

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser("Build predition model from data")
    parser.add_argument("-data_annotation")
    parser.add_argument("-variant_annotation")
    parser.add_argument("-model_input")
    parser.add_argument("-output")
    parser.add_argument("-parsimony", type=int, default=logging.INFO)
    args = parser.parse_args()

    Logging.configure_logging(args.parsimony)
    run(args)