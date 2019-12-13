__author__ = "alvaro barbeira"

import subprocess
import re
import os
import logging
import gzip
import numpy
import pandas
import collections
import shutil
import traceback

import pyarrow.parquet as pq

from genomic_tools_lib import Utilities, Logging
from genomic_tools_lib.data_management import TextFileTools
from genomic_tools_lib.individual_data import Utilities as StudyUtilities
from genomic_tools_lib.file_formats import Parquet
from genomic_tools_lib.miscellaneous import matrices, Genomics

def _intermediate_folder(intermediate_folder, gene): return os.path.join(intermediate_folder, gene)
def _y_folder(intermediate_folder, gene): return os.path.join(_intermediate_folder(intermediate_folder, gene), "y")
def _x_path(intermediate_folder, gene): return os.path.join(_intermediate_folder(intermediate_folder, gene), "x.txt")
def _info_path(intermediate_folder, gene): return os.path.join(_intermediate_folder(intermediate_folder, gene), "info.txt")
def _outdir(intermediate_folder, gene): return os.path.join(_intermediate_folder(intermediate_folder, gene), "out")
def _weights(intermediate_folder, gene): return os.path.join(_outdir(intermediate_folder, gene), "result.m.est")
def _summary(intermediate_folder, gene): return os.path.join(_outdir(intermediate_folder, gene), "result.m.stats")
def _execution_script(intermediate_folder, gene): return os.path.join(_intermediate_folder(intermediate_folder, gene), "ctimp.sh")

########################################################################################################################
def prepare_ctimp(script_path, seed, intermediate_folder, data_annotation_, features_, features_data_, d_):
    _i = _intermediate_folder(intermediate_folder, data_annotation_.gene_id)
    if os.path.exists(_i):
        logging.info("intermediate folder for %s exists, aborting", data_annotation_.gene_id)
        raise RuntimeError("Dirty folder")
    os.makedirs(_i)

    save_expression(intermediate_folder, data_annotation_.gene_id, d_, features_data_)
    save_x(intermediate_folder, data_annotation_.gene_id, features_, features_data_)
    execution_script(script_path, seed, intermediate_folder, data_annotation_.gene_id)

def save_x(intermediate_folder, gene, features_, features_data_):
    Utilities.save_dataframe(
        features_data_.drop("individual", axis=1),
        _x_path(intermediate_folder, gene),
        header=False, sep=" ")

    Utilities.save_dataframe(
        features_[["id", "allele_0", "allele_1"]].rename(columns={"id":"SNP", "allele_0":"REF.0.", "allele_1":"ALT.1."}),
        _info_path(intermediate_folder, gene))

def save_expression(intermediate_folder, gene, d_, features_data_):
    y_folder = _y_folder(intermediate_folder, gene)
    os.makedirs(_y_folder(intermediate_folder, gene))
    for k,v in d_.items():
        if not gene in v:
            logging.log(8, "%s not present in %s", gene, k)
            continue
        p = os.path.join(y_folder,k) + ".txt"
        v = v.merge(features_data_[["individual", "id"]], on="individual")[["id", gene]]
        Utilities.save_dataframe(v, p, header=False)

def execution_script(script_path, seed, intermediate_folder, gene):
    r_ = os.path.split(script_path)[0]
    script=\
"""#!/bin/bash

cd {run_path}

#mkdir -p {outdir}

Rscript {script_path} \\
-dosage {dosage} \\
-info {info} \\
-expression {expression} \\
-ntune 50 \\
-nfold 5 \\
-outdir {outdir} \\
-seed {seed} \\
-gene_id result #> /dev/null
""".format(run_path=r_, script_path=script_path, seed=seed,
        dosage=_x_path(intermediate_folder, gene), info=_info_path(intermediate_folder, gene),
        expression=_y_folder(intermediate_folder, gene), outdir=_outdir(intermediate_folder, gene))
    with open(_execution_script(intermediate_folder,gene), "w") as s:
        s.write(script)

def setup_output(output_prefix, tissue_names, WEIGHTS_FIELDS, SUMMARY_FIELDS):
    weights = {}
    summaries = {}
    covariances = {}
    for t in tissue_names:
        w_ = "{}_{}_t_weights.txt.gz".format(output_prefix, t)
        if os.path.exists(w_):
            raise RuntimeError("weights exist! delete them or move them")
        weights[t] = gzip.open(w_, "w")
        weights[t].write(("\t".join(WEIGHTS_FIELDS) + "\n").encode())

        summaries[t] = gzip.open("{}_{}_t_summary.txt.gz".format(output_prefix, t), "w")
        summaries[t].write(("\t".join(SUMMARY_FIELDS) + "\n").encode())

        covariances[t] = gzip.open("{}_{}_t_covariance.txt.gz".format(output_prefix, t), "w")
        covariances[t].write("GENE RSID1 RSID2 VALUE\n".encode())
    return weights, summaries, covariances

def set_down(weights, summaries, covariances, tissue_names, failed_run):
    for t in tissue_names:
        weights[t].close()
        summaries[t].close()
        covariances[t].close()
        if failed_run:
            os.remove(os.path.realpath(weights[t].name))
            os.remove(os.path.realpath(summaries[t].name))
            os.remove(os.path.realpath(covariances[t].name))

########################################################################################################################
def run(args):
    Utilities.maybe_create_folder(args.intermediate_folder)
    Utilities.ensure_requisite_folders(args.output_prefix)

    logging.info("Opening data")
    p_ = re.compile(args.data_name_pattern)
    f = [x for x in sorted(os.listdir(args.data_folder)) if p_.search(x)]
    tissue_names = [p_.search(x).group(1) for x in f]
    data = []
    for i in range(0,len(tissue_names)):
        logging.info("Loading %s", tissue_names[i])
        data.append((tissue_names[i], pq.ParquetFile(os.path.join(args.data_folder, f[i]))))
    data = collections.OrderedDict(data)
    available_data = {x for p in data.values() for x in p.metadata.schema.names}

    logging.info("Preparing output")
    WEIGHTS_FIELDS=["gene", "rsid", "varID", "ref_allele", "eff_allele", "weight"]
    SUMMARY_FIELDS=["gene", "genename", "gene_type", "alpha", "n_snps_in_window", "n.snps.in.model", "rho_avg", "pred.perf.R2", "pred.perf.pval"]

    Utilities.ensure_requisite_folders(args.output_prefix)

    if args.skip_regression:
        weights, summaries, covariances = None, None, None
    else:
        weights, summaries, covariances = setup_output(args.output_prefix, tissue_names, WEIGHTS_FIELDS, SUMMARY_FIELDS)

    logging.info("Loading data annotation")
    data_annotation = StudyUtilities._load_gene_annotation(args.data_annotation)
    data_annotation = data_annotation[data_annotation.gene_id.isin(available_data)]
    if args.chromosome or (args.sub_batches and args.sub_batch):
        data_annotation = StudyUtilities._filter_gene_annotation(data_annotation, args.chromosome, args.sub_batches, args.sub_batch)
    logging.info("Kept %i entries", data_annotation.shape[0])

    logging.info("Opening features annotation")
    if not args.chromosome:
        features_metadata = pq.read_table(args.features_annotation).to_pandas()
    else:
        features_metadata = pq.ParquetFile(args.features_annotation).read_row_group(args.chromosome-1).to_pandas()

    if args.discard_palindromic_snps:
        logging.info("Discarding palindromic snps")
        features_metadata = Genomics.discard_gtex_palindromic_variants(features_metadata)

    if args.chromosome and args.sub_batches:
        logging.info("Trimming variants")
        features_metadata = StudyUtilities.trim_variant_metadata_on_gene_annotation(features_metadata, data_annotation, args.window)

    if args.rsid_whitelist:
        logging.info("Filtering features annotation")
        whitelist = TextFileTools.load_list(args.rsid_whitelist)
        whitelist = set(whitelist)
        features_metadata = features_metadata[features_metadata.rsid.isin(whitelist)]

    logging.info("Opening features")
    features = pq.ParquetFile(args.features)

    logging.info("Setting R seed")
    seed = numpy.random.randint(1e8)

    if args.run_tag:
        d = pandas.DataFrame({"run":[args.run_tag], "cv_seed":[seed]})[["run", "cv_seed"]]
        for t in tissue_names:
            Utilities.save_dataframe(d, "{}_{}_t_runs.txt.gz".format(args.output_prefix, t))

    failed_run=False
    try:
        for i, data_annotation_ in enumerate(data_annotation.itertuples()):
            logging.log(9, "processing %i/%i:%s", i+1, data_annotation.shape[0], data_annotation_.gene_id)
            logging.log(8, "loading data")
            d_ = {}
            for k, v in data.items():
                d_[k] = Parquet._read(v, [data_annotation_.gene_id], to_pandas=True)
            features_ = Genomics.entries_for_gene_annotation(data_annotation_, args.window,
                                                             features_metadata)

            if features_.shape[0] == 0:
                logging.log(9, "No features available")
                continue

            features_data_ = Parquet._read(features, [x for x in features_.id.values], to_pandas=True)
            features_data_["id"] = range(1, features_data_.shape[0] + 1)
            features_data_ = features_data_[["individual", "id"] + [x for x in features_.id.values]]

            logging.log(8, "training")
            prepare_ctimp(args.script_path, seed, args.intermediate_folder, data_annotation_, features_, features_data_, d_)
            del(features_data_)
            del(d_)
            if args.skip_regression:
                continue

            subprocess.call(["bash", _execution_script(args.intermediate_folder, data_annotation_.gene_id)])

            w = pandas.read_table(_weights(args.intermediate_folder, data_annotation_.gene_id), sep="\s+")
            s = pandas.read_table(_summary(args.intermediate_folder, data_annotation_.gene_id), sep="\s+")

            for e_, entry in enumerate(s.itertuples()):
                entry_weights = w[["SNP", "REF.0.", "ALT.1.", entry.tissue]].rename(columns={
                    "SNP":"varID", "REF.0.":"ref_allele", "ALT.1.":"eff_allele", entry.tissue:"weight"
                })
                entry_weights = entry_weights[entry_weights.weight != 0]
                entry_weights = entry_weights.assign(gene = data_annotation_.gene_id)
                entry_weights = entry_weights.merge(features_, left_on="varID", right_on="id", how="left")
                entry_weights = entry_weights[WEIGHTS_FIELDS]
                if args.output_rsids:
                    entry_weights.loc[entry_weights.rsid == "NA", "rsid"] = entry_weights.loc[entry_weights.rsid == "NA", "varID"]
                weights[entry.tissue].write(entry_weights.to_csv(sep="\t", index=False, header=False, na_rep="NA").encode())

                entry_summary = s[s.tissue == entry.tissue].rename(columns={"zscore_pval":"pred.perf.pval", "rho_avg_squared":"pred.perf.R2"})
                entry_summary = entry_summary.assign(gene = data_annotation_.gene_id, alpha=0.5,
                    genename = data_annotation_.gene_name, gene_type= data_annotation_.gene_type, n_snps_in_window = features_.shape[0])
                entry_summary["n.snps.in.model"] = entry_weights.shape[0]
                #must repeat strings beause of weird pandas indexing issue
                entry_summary = entry_summary.drop(["R2", "n", "tissue"], axis=1)[["gene", "genename", "gene_type", "alpha", "n_snps_in_window", "n.snps.in.model", "rho_avg", "pred.perf.R2", "pred.perf.pval"]]
                summaries[entry.tissue].write(entry_summary.to_csv(sep="\t", index=False, header=False, na_rep="NA").encode())

                features_data_ = Parquet._read(features, [x for x in entry_weights.varID.values], to_pandas=True)
                var_ids = [x for x in entry_weights.varID.values]
                cov = numpy.cov([features_data_[k] for k in var_ids], ddof=1)
                ids = [x for x in entry_weights.rsid.values] if args.output_rsids else var_ids
                cov = matrices._flatten_matrix_data([(data_annotation_.gene_id, ids, cov)])
                for cov_ in cov:
                    l = "{} {} {} {}\n".format(cov_[0], cov_[1], cov_[2], cov_[3]).encode()
                    covariances[entry.tissue].write(l)

            if not args.keep_intermediate_folder:
                logging.info("Cleaning up")
                shutil.rmtree(_intermediate_folder(args.intermediate_folder, data_annotation_.gene_id))

            if args.MAX_M and i >= args.MAX_M:
                logging.info("Early abort")
                break

    except Exception as e:
        logging.info("Exception running model training:\n%s", traceback.format_exc())
        failed_run=True
    finally:
        pass
        # This is wrong but I must thing about what to do
        #if not args.keep_intermediate_folder:
        #    shutil.rmtree(args.intermediate_folder)

    if not args.skip_regression:
        set_down(weights, summaries, covariances, tissue_names, failed_run)
    if failed_run:
        logging.info("Errors found")
    else:
        logging.info("Finished")


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser("Train Elastic Net prediction models from GLMNET")
    parser.add_argument("-run_tag")
    parser.add_argument("-script_path")
    parser.add_argument("-features")
    parser.add_argument("-features_annotation")
    parser.add_argument("-data_folder")
    parser.add_argument("-data_name_pattern")
    parser.add_argument("-data_annotation")
    parser.add_argument("-window", type = int)
    parser.add_argument("-intermediate_folder")
    parser.add_argument("--output_rsids", action="store_true")
    parser.add_argument("--chromosome", type = int)
    parser.add_argument("--sub_batches", type = int)
    parser.add_argument("--sub_batch", type =int)
    parser.add_argument("--rsid_whitelist")
    parser.add_argument("--keep_intermediate_folder", action="store_true")
    parser.add_argument("--MAX_M", type=int)
    parser.add_argument("--skip_regression", action="store_true")
    parser.add_argument("--discard_palindromic_snps", action="store_true")
    parser.add_argument("-output_prefix")
    parser.add_argument("-parsimony", default=10, type=int)
    args = parser.parse_args()
    Logging.configure_logging(args.parsimony)

    run(args)