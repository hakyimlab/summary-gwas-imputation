#!/usr/bin/env python
__author__ = "alvaro barbeira"

import logging
from timeit import default_timer as timer

import numpy

from genomic_tools_lib import Logging, Utilities, Exceptions
from genomic_tools_lib.data_management import DataFrameStreamer, GTExMisc
from genomic_tools_lib.file_formats import Parquet
from genomic_tools_lib.file_formats.eqtl import GTEx

def _skip_eqtl(comps):
    """Deprecated."""
    #skip anything not biallelic
    variant = comps[1]
    c = variant.split("_")
    return len(c[2]) != 1 or len(c[3]) != 1

def _skip_missing_in_key(comps, key_to_snp):
    """Deprecated"""
    return not comps[1] in key_to_snp

def _process(d, key_to_snp, how="left"):
    k = [(k_, key_to_snp[k_]) for k_ in d.variant_id if k_ in key_to_snp]
    k = Utilities.to_dataframe(k, ["variant_id", "rsid"])
    d = d.merge(k, on="variant_id", how=how)
    d = d.rename(columns={"gene_id":"gene", "pval_nominal":"pvalue", "slope":"beta", "slope_se":"se"})
    d = d[list(GTEx.GTExAllAssociations._fields)]
    d = d.assign(maf = d.maf.astype(numpy.float32), beta = d.beta.astype(numpy.float32), se = d.se.astype(numpy.float32))
    return d

def load_annotation(args):
    if args.snp_annotation:
        logging.info("Loading SNP annotation file for variant-rsid mapping")
        return GTExMisc.load_gtex_variant_to_rsid(args.snp_annotation)
    elif args.snp_annotation_from_parquet_metadata:
        logging.info("Loading Parquet metadata for variant-rsid mapping")
        return Parquet.variant_key_value_from_metadata(args.snp_annotation_from_parquet_metadata)

    raise Exceptions.ReportableException("Provide a file annotation")

def _data_sink(args):
    return Parquet.ParquetDataFrameSink(args.parquet_output, GTEx.pyarrow_schema)

def run(args):
    start = timer()
    Utilities.ensure_requisite_folders(args.parquet_output)

    logging.info("Loading snp annotation")
    key_to_snp = load_annotation(args)

    logging.info("Processing eqtl")
    _skip = (lambda x: _skip_missing_in_key(x, key_to_snp)) if args.restrict_to_annotation else None
    streamer_ = DataFrameStreamer.data_frame_streamer(args.gtex_eqtl_file, sanitize=True,
                    to_numeric=["maf", "pval_nominal", "slope", "slope_se"], sentinel_column="gene_id", additional_skip_row_check=_skip)

    with _data_sink(args) as sink:
        for i,d in enumerate(streamer_):
            if d.shape[0] == 0:
                logging.log(8, "Skipping %d", i)
                continue
            logging.log(8, "Processing %d/%s", i, d.gene_id[0])
            p = _process(d, key_to_snp)
            sink.sink(p)

    end = timer()
    logging.info("Ran conversion in %s", str(end-start))

if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser("Convert GTEx association file")
    parser.add_argument("-snp_annotation", help="File describing snps")
    parser.add_argument("-gtex_eqtl_file", help="Which GTEx file to convert")
    parser.add_argument("-snp_annotation_from_parquet_metadata", help="Load a genotype study metadata (parquet format) for the required information")
    parser.add_argument("-restrict_to_annotation", help="keep only snps in the annotation", action="store_true")
    parser.add_argument("-parquet_output", help="Where to save the output")
    parser.add_argument("-parsimony", help="Log parsimony level. 1 is everything being logged. 10 is only high level messages, above 10 will hardly log anything", default = 10, type=int)

    args = parser.parse_args()

    Logging.configure_logging(args.parsimony)

    run(args)