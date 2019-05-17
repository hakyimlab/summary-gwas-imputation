__author__ = "alvaro barbeira"
import logging
import os
import re
import numpy
import pandas
import traceback

from genomic_tools_lib import Utilities, Logging
from genomic_tools_lib.data_management import DataFrameStreamer, TextFileTools
from genomic_tools_lib.file_formats.Parquet import ParquetSingleSplitStudy

import rpy2.robjects as robjects
import rpy2.robjects.numpy2ri
import rpy2.robjects.pandas2ri
import rpy2.rinterface
rpy2.robjects.numpy2ri.activate()


def get_study(parquet_folder, parquet_pattern, parquet_metadata):
    import pyarrow as pa
    import pyarrow.parquet as pq

    logging.info("Reading variant metadata")
    #vm = pq.ParquetFile(parquet_metadata).read().to_pandas()
    v = pq.ParquetFile(parquet_metadata).read(["id"])
    v = set(v.column(0).to_pylist())

    logging.info("Acquiring genotype files")
    r = re.compile(parquet_pattern)
    names = os.listdir(args.parquet_genotype_folder)
    files = {r.search(x).group(1):os.path.join(parquet_folder, x) for x in names if r.search(x)}
    #files = {r.search(x).group(1):pq.ParquetFile(os.path.join(parquet_folder,x)) for x in names if r.search(x)}

    study = ParquetSingleSplitStudy(files, None)
    return study, v

def _dump(p, d, cov):
    Utilities.save_dataframe(d, p+ "_d.txt.gz")
    import gzip
    with gzip.open(p+"_m.txt.gz", "w") as f:
        for i in cov:
            f.write("{}\n".format("\t".join(map(str,i))).encode())

from rpy2.robjects.packages import importr
susieR = importr('susieR')
susie_bhat = susieR.susie_bhat
susie_z = susieR.susie_z
summary = robjects.r["summary"]

def _do_susie(d, study, variants_whitelist, n, specific_individuals, mode):
    d_ = d[d.variant_id.isin(variants_whitelist)]
    variants = [x for x in d_.variant_id.values]
    if specific_individuals:
        X = study.get_variants(variants, specific_individuals=specific_individuals)
        del X["individual"]
    else:
        X = study.get_variants(variants, omit_individuals=True)
    # present_variants = set(X.keys())
    # d_ = d_[d_.variant_id.isin(present_variants)]
    X = [X[x] for x in d_.variant_id]
    cov = numpy.cov(X, ddof=1)
    cov = robjects.r.matrix(robjects.FloatVector(cov.flatten()), nrow=len(variants))

    r = None
    if mode is None or mode=="bhat":
        b_ = robjects.FloatVector(d_.slope.values)
        s_ = robjects.FloatVector(d_.slope_se.values)
        r = susie_bhat(b_, s_, cov, n=n)
    elif mode == "z":
        z_ = robjects.FloatVector(d_.slope/d_.slope_se)
        r = susie_z(z_, cov)
    else:
        raise RuntimeError("Unrecognized mode")
    return r,d_

def _void_cs(status=None):
    return pandas.DataFrame({"cs": [None], "cs_log10bf": [None], "cs_avg_r2": [None], "cs_min_r2": [None], "variable": [None], "status": [status], "var_id": [None]})

def _void_var():
    return pandas.DataFrame({"variable":[], "variable_prob":[], "cs":[]})

def _process_result(res_, d_, gene):
    sum_ = summary(res_)
    if sum_[1] is rpy2.rinterface.NULL:
        logging.log(9, "No credible set")
        cs = _void_cs("no_credible_set")
    else:
        cs = rpy2.robjects.pandas2ri.ri2py(sum_[1]).assign(status=None)
        csids = []
        for t in cs.itertuples():
            s = map(int, t.variable.split(","))
            s = map(lambda x: d_.iloc[x - 1].variant_id, s)
            csids.append(",".join(s))

    vars = rpy2.robjects.pandas2ri.ri2py(sum_[0])
    pp = numpy.sum(vars.variable_prob)
    vars = vars[vars.cs > -1]
    var_ids = [d_.iloc[int(x) - 1].variant_id for x in vars.variable]

    cs = cs.assign(gene_id=gene, pp_sum=pp)
    vars = vars.assign(gene_id=gene, var_id=var_ids)
    return cs, vars

def run(args):
    if os.path.exists(args.cs_output) or os.path.exists(args.var_output):
        logging.info("Output exists. Nope.")
        return

    study, variants_whitelist = get_study(args.parquet_genotype_folder, args.parquet_genotype_pattern, args.parquet_genotype_metadata)

    #_skip = lambda x: x not in variants_whitelist
    columns = ["maf", "pval_nominal", "slope", "slope_se"]
    eqtl_streamer = DataFrameStreamer.data_frame_streamer(args.eqtl, sanitize=True, to_numeric=columns, sentinel_column="gene_id")

    individuals = None if not args.restrict_to_individuals else TextFileTools.load_list(args.restrict_to_individuals)

    genes = None if not args.restrict_to_genes else set(TextFileTools.load_list(args.restrict_to_genes))

    cs_results = []
    var_results = []
    logging.info("Beggining process")
    MAX_N=args.MAX_N
    n=args.sample_size
    for i, d in enumerate(eqtl_streamer):
        if MAX_N and i > MAX_N:
            logging.info("Early exit")
            break
        gene = d.gene_id.values[0]
        if genes is not None and gene.split('.')[0] not in genes:
            logging.log(9, "Skipping gene: %s", gene)
            continue
        logging.log(9, "Processing gene %i:%s", i+1, gene)
        d = d.loc[(~d.slope_se.isnull()) & (d.slope!=0) & (~d.slope.isnull())]
        try:
            res_, d_ = _do_susie(d, study, variants_whitelist, n, individuals, args.mode)
            cs, vars =_process_result(res_, d_, gene)
        except Exception as e:
            logging.log(9, "Error while doing susie:\n%s", traceback.format_exc())
            cs = _void_cs("susie_error").assign(gene_id=gene, pp_sum=None)
            vars = _void_var().assign(gene_id=[gene], var_id=[None])

        cs_results.append(cs)
        #if vars.shape[1]>0:
        var_results.append(vars)

    if len(cs_results) > 0:
        logging.info("Saving")
        cs_results = pandas.concat(cs_results)[["gene_id", "cs", "cs_avg_r2", "cs_log10bf", "cs_min_r2", "var_id", "pp_sum", "status"]]
        Utilities.ensure_requisite_folders(args.cs_output)
        Utilities.save_dataframe(cs_results, args.cs_output)
    else:
        logging.info('No results')

    if len(var_results) > 0:
        var_results = pandas.concat(var_results)[["gene_id", "var_id",  "cs", "variable_prob"]]
        Utilities.ensure_requisite_folders(args.var_output)
        Utilities.save_dataframe(var_results, args.var_output)
    logging.info("Ran susie")

if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser("Run COLOC")
    parser.add_argument("-parquet_genotype_folder", help="Parquet Genotype file")
    parser.add_argument("-parquet_genotype_pattern", help="Parquet Genotype file")
    parser.add_argument("-parquet_genotype_metadata", help="Parquet Genotype variant metadata file")
    parser.add_argument("-restrict_to_individuals", help="filter to individuals")
    parser.add_argument("-restrict_to_genes", help="filter to genes")
    parser.add_argument("--mode", help="'bhat' or 'z' (bhat is default)")
    parser.add_argument("-eqtl", help="Run on a GTEX-like eqtl summary stats file")
    parser.add_argument("-sample_size", help="number of samples", type=int)
    parser.add_argument("-cs_output", help="Credible sets.")
    parser.add_argument("-var_output", help="variables")
    parser.add_argument("-parsimony", help="Log verbosity level. 1 is everything being logged. 10 is only high level messages, above 10 will hardly log anything", default = "10")
    parser.add_argument("-MAX_N", type=int)

    args = parser.parse_args()

    Logging.configure_logging(int(args.parsimony), with_date=True)

    run(args)
