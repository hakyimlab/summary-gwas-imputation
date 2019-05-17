__author__ = "alvaro barbeira"

import os
import logging
import pandas
import sqlite3

import pyarrow.parquet as pq

from genomic_tools_lib import Utilities, Logging
from genomic_tools_lib.individual_data import Utilities as StudyUtilities
from genomic_tools_lib.miscellaneous import Models
from genomic_tools_lib.file_formats import Miscellaneous

def get_weights(x_weights, id_whitelist=None):
    if x_weights[1] == "PIP":
        w = Miscellaneous.dapg_signals(x_weights[0], float(x_weights[2]), id_whitelist)
        w = w.rename(columns={"gene":"gene_id", "pip":"w", "variant_id":"id"})
    else:
        raise RuntimeError("unsupported weights argument")
    return w

def run(args):
    if os.path.exists(args.output):
        logging.info("output exists already, delete it or move it")
        return

    logging.info("Starting")
    Utilities.ensure_requisite_folders(args.output)

    logging.info("Loading data annotation")
    gene_annotation = StudyUtilities.load_gene_annotation(args.gene_annotation)
    gene_annotation = gene_annotation.rename({"gene_name":"genename"}, axis=1)[["gene_id", "genename", "gene_type"]]

    logging.info("Loading variant annotation")
    features_metadata = pq.read_table(args.features_annotation).to_pandas()

    logging.info("Loading spec")
    weights = get_weights(args.spec)

    w = weights.merge(features_metadata[["id", "allele_0", "allele_1", "rsid"]], on="id", how="left")
    w = w.rename({"allele_0":"ref_allele", "allele_1":"eff_allele", "id":"varID"}, axis=1)
    w["gene"] = w.gene_id.str.cat(w.cluster_id.astype(str), sep="_")
    w = w.drop(["w", "cluster_id"], axis=1)
    w = w.sort_values(by="gene").assign(weight = 1)

    logging.info("Building models")
    with sqlite3.connect(args.output) as conn:
        w.drop("gene_id", axis=1).fillna("NA")[["gene", "rsid", "varID", "ref_allele", "eff_allele", "weight"]].to_sql("weights", conn, index=False)

        e = w[["gene_id", "gene"]].merge(gene_annotation, on="gene_id").drop("gene_id", axis=1)
        e["n_snps_in_window"] = None
        e["n.snps.in.model"] = 1
        e["pred.perf.pval"] = None
        e["pred.perf.qval"] = None
        e["pred.perf.R2"] = None
        e = e[["gene", "genename", "gene_type", "n_snps_in_window", "n.snps.in.model", "pred.perf.R2", "pred.perf.pval", "pred.perf.qval"]]

        e.to_sql("extra", conn, index=False)

        Models.model_indexes(conn)

    logging.info("Finished")

if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser("Train Elastic Net prediction models from GLMNET")
    parser.add_argument("-spec", nargs="+")
    parser.add_argument("-gene_annotation")
    parser.add_argument("-features_annotation")
    parser.add_argument("-output")
    parser.add_argument("-parsimony", default=10, type=int)
    args = parser.parse_args()
    Logging.configure_logging(args.parsimony)

    run(args)