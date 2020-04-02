#!/usr/bin/env python
__author__ = "alvaro barbeira"
import os
import logging
import numpy
import pandas
import sqlite3
import gzip

from genomic_tools_lib import Utilities, Logging
from genomic_tools_lib.file_formats import Gencode
from genomic_tools_lib.miscellaneous import Models
from genomic_tools_lib.data_management import DataFrameStreamer, KeyedDataSource

def _gene_annotation(path):
    genes = Gencode.load(path, selected_key_value_pairs = [Gencode.GFTF.K_GENE_ID, Gencode.GFTF.K_GENE_NAME, Gencode.GFTF.K_GENE_TYPE])
    genes_ = {}
    types_ = {}
    for t in genes.itertuples():
        genes_[t.gene_id] = t.gene_name
        types_[t.gene_id] = t.gene_type
    return genes_, types_

def run(args):
    if os.path.exists(args.output):
        logging.info("Output already exists, either delete it or move it")
        return

    logging.info("Loading snp names")
    snps = KeyedDataSource.load_data(args.snp_annotation, "variant_id", args.rsid_column)

    logging.info("Loading gene annotation")
    genes, types = _gene_annotation(args.gene_annotation)

    with sqlite3.connect(args.output) as conn:
        logging.info("Processing")

        streamer = DataFrameStreamer.data_frame_streamer(args.input,
            header=["tissue_name", "gene_id", "variant_id", "weight", "beta", "se"],
            to_numeric=["weight", "beta", "se"], sentinel_column="gene_id")
        extra = []
        for i, d in enumerate(streamer):
            g_ = d.gene_id.values[0]
            logging.log(9, "processing %i:%s", i+1, g_)
            d = d.loc[d.weight != 0]
            if args.snp_zscore_threshold:
                d = d.assign(zscore=numpy.abs(d.beta / d.se))
                d = d.loc[d.zscore > args.snp_zscore_threshold]

            if d.shape[0] == 0:
                logging.log(9, "no good snps left")
                continue

            extra.append((g_, genes[g_], types[g_], d.shape[0], numpy.nan, numpy.nan, numpy.nan))

            d = d[["gene_id", "variant_id", "weight"]].rename(columns={"gene_id":"gene", "variant_id":"varID"})
            effect, non_effect, rsid = [], [], []
            for t in d.itertuples():
                c_ = t.varID.split("_")
                effect.append(c_[3])
                non_effect.append(c_[2])
                r_ = snps[t.varID]
                rsid.append(r_ if r_ != "." else t.varID)
            d = d.assign(ref_allele = non_effect, eff_allele = effect, rsid = rsid)[["gene", "rsid", "varID", "ref_allele", "eff_allele", "weight"]]
            d.to_sql("weights", conn, index=False, if_exists="append")

        extra = pandas.DataFrame(extra, columns=["gene", "genename", "gene_type", "n.snps.in.model", "pred.perf.R2","pred.perf.pval", "pred.perf.qval"])
        extra.to_sql("extra", conn, index=False)

        logging.info("Creating indices")
        Models.model_indexes(conn)

    logging.info("Finished building model.")


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser("Convert -twas builder- output to Predixcan models")

    parser.add_argument("-input")
    parser.add_argument("-output")
    parser.add_argument("-snp_annotation", help="gtex-like snp annotation")
    parser.add_argument("-gene_annotation", help="gtf-like or gencode file")
    parser.add_argument("-snp_zscore_threshold", help="Optional. Keep only snps with good enough (absolute value) marginal zscore", default=None, type=float)
    parser.add_argument("-parsimony", help="Log verbosity level. 1 is everything being logged. 10 is only high level messages, above 10 will hardly log anything", default=10, type=int)
    parser.add_argument("-rsid_column", help = "Column name with rsid variant identifiers", default = "rs_id_dbSNP150_GRCh38p7")
    args =parser.parse_args()
    Logging.configure_logging(args.parsimony)
    run(args)
