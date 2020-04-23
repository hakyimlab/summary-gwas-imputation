__author__ = "alvaro barbeira"

import sys
import os
import logging
import gzip
import collections
import re

import numpy
import pandas
import scipy
import math
import statsmodels.formula.api as smf
from sklearn.model_selection import KFold
from sklearn.metrics import r2_score

import pyarrow.parquet as pq

from genomic_tools_lib import Utilities, Logging
from genomic_tools_lib.data_management import TextFileTools
from genomic_tools_lib.individual_data import Utilities as StudyUtilities
from genomic_tools_lib.file_formats import Parquet, Miscellaneous
from genomic_tools_lib.miscellaneous import matrices, Genomics, Math

###############################################################################
import rpy2.robjects as robjects
from rpy2.robjects import numpy2ri
numpy2ri.activate()

from rpy2.robjects import pandas2ri
pandas2ri.activate()

def initialize():
    global train_elastic_net
    global set_seed
    path = os.path.split(sys.argv[0])[0]
    path = os.path.join(path, "elastic_net.R")
    robjects.r['source'](path)
    train_elastic_net = robjects.r['train_elastic_net']
    set_seed = robjects.r['set_seed']
###############################################################################

def _save(d_, features_, features_data_, gene):
    pandas.DataFrame({gene:d_[gene]}).to_csv("y.txt", index=False, sep="\t")
    pandas.DataFrame(collections.OrderedDict([(v,features_data_[v]) for v in features_.id.values])).to_csv("x.txt", index=False, sep="\t")

def get_weights(x_weights, id_whitelist):
    if x_weights[1] == "PIP":
        w = Miscellaneous.dapg_signals(x_weights[0], float(x_weights[2]), id_whitelist)
        w = w.rename(columns={"gene":"gene_id", "pip":"w", "variant_id":"id"})
        w.w = 1 - w.w #Less penalty to the more probable snps
    else:
        raise RuntimeError("unsupported weights argument")
    return w

def get_dapg_preparsed(weights):
    df = Miscellaneous.dapg_preparsed(weights)
    df.w = 1 - df.w
    return df

###############################################################################

def train_elastic_net_wrapper(features_data_, features_, d_, data_annotation_, x_w=None, prune=True, nested_folds=10):
    x = numpy.array([features_data_[v] for v in features_.id.values])
    dimnames = robjects.ListVector(
        [(1, robjects.StrVector(d_["individual"])), (2, robjects.StrVector(features_.id.values))])
    x = robjects.r["matrix"](robjects.FloatVector(x.flatten()), ncol=features_.shape[0], dimnames=dimnames)
    y = robjects.FloatVector(d_[data_annotation_.gene_id])
    nested_folds = robjects.FloatVector([nested_folds])
    #py2ri chokes on None.
    if x_w is None:
        res = train_elastic_net(y, x, n_train_test_folds=nested_folds)
    else:
        res = train_elastic_net(y, x, penalty_factor=x_w, n_train_test_folds=nested_folds)  # observation weights, not explanatory variable weight :( , x_weight = x_w)
    return pandas2ri.ri2py(res[0]), pandas2ri.ri2py(res[1])

###############################################################################

def ols_pred_perf(data, n_folds=10):
    kf = KFold(n_splits=n_folds, shuffle=True)

    rho_f=[]
    R2_f=[]
    zscore_f=[]
    for train_index, test_index in kf.split(data):
        train_ = data.iloc[train_index]
        fold_= smf.ols('y ~ {}'.format(" + ".join([x for x in train_.columns if x !="y"])), data=train_).fit()
        test_ = data.iloc[test_index]
        y_predicted = fold_.predict(test_)
        if numpy.std(y_predicted) != 0:
            score = r2_score(test_.y, y_predicted)
            rho = numpy.corrcoef(test_.y, y_predicted)[0,1]
            zscore = numpy.arctanh(rho)*numpy.sqrt(len(y_predicted) - 3)
        else:
            score = 0
            rho = 0
            zscore = 0
        R2_f.append(score)
        rho_f.append(rho)
        zscore_f.append(zscore)

    rho_avg = numpy.average(rho_f)
    zscore_est = numpy.sum(zscore_f)/numpy.sqrt(n_folds)
    zscore_pval = scipy.stats.norm.sf(zscore_est)
    d = {"test_R2_avg": [numpy.average(R2_f)], "test_R2_sd": [numpy.std(R2_f)],
         "rho_avg": [rho_avg], "rho_avg_squared": [rho_avg**2], "rho_se":[numpy.std(rho_f)],
         "rho_zscore":[zscore_est], "zscore_pval": [zscore_pval], "nested_cv_fisher_pval":[None], "nested_cv_converged":n_folds}
    return pandas.DataFrame(d)

def prune(data):
    if data.shape[1] == 1:
        return data
    cor = numpy.corrcoef(data.values.T)
    discard=set()
    for i in range(0, cor.shape[0]):
        for j in range(i, cor.shape[1]):
            if i==j:
                continue
            if i in discard:
                continue
            if math.abs(cor[i][j]) >= 0.95:
                discard.add(j)
    discard = data.columns[list(discard)].values
    return data.drop(discard, axis=1)


def train_ols(features_data_, features_, d_, data_annotation_, x_w=None, prune=True, nested_folds=10):
    ids=[]
    data = {}
    for v in features_.id.values:
        x = Math.standardize(features_data_[v])
        if x is not None:
            data[v] = x
            ids.append(v)
    data = pandas.DataFrame(data)

    if prune:
        data = prune(data)
    ids = data.columns.values
    if len(ids) == 0:
        w = pandas.DataFrame({"feature":[], "weight":[]})
        s = pandas.DataFrame({"test_R2_avg": [], "test_R2_sd": [],
         "rho_avg": [], "rho_avg_squared": [], "rho_se":[],
         "rho_zscore":[], "zscore_pval": [], "nested_cv_fisher_pval":[],
        "alpha":[], "n_snps_in_window":[], "cv_R2_avg":[], "cv_R2_sd":[], "in_sample_R2":[], "n.snps.in.model":[]})
        return w,s

    data["y"] = Math.standardize(d_[data_annotation_.gene_id])

    results = smf.ols('y ~ {}'.format(" + ".join(ids)), data=data).fit()
    weights = results.params[1:].to_frame().reset_index().rename(columns={"index": "feature", 0: "weight"})
    summary = ols_pred_perf(data, nested_folds)
    summary = summary.assign(alpha=None, n_snps_in_window=features_.shape[0],
                           cv_R2_avg=None, cv_R2_sd=None, in_sample_R2=None)
    summary["n.snps.in.model"] = len(ids)
    return weights, summary


########################################################################################################################

def process(w, s, c, data_handler, data_annotation_, features_handler,
             weights, summary_fields, train, postfix=None,
            nested_folds=10):
    """

    :param w: weights file handle
    :param s: summary file handle
    :param c: covariance file handle
    :param data_handler: class DataHandler
    :param data_annotation_: Pandas tuple with attributes 'gene_name', 'gene_id', 'gene_type'
    :param features_handler: class FeaturesHandler
    :param weights: Bool
    :param summary_fields: list of strings
    :param train: function. Training function
    :param postfix: int. If doing repeats, this is the repeat number.
    :param nested_folds: int. Number of nested folds.
    :return:
    """
    gene_id_ = data_annotation_.gene_id if postfix is None else "{}-{}".format(data_annotation_.gene_id, postfix)
    logging.log(8, "Loading phenotype data")
    d_ = data_handler.load_pheno(data_annotation_.gene_id)
    # d_ = Parquet._read(data, [data_annotation_.gene_id])

    features_ = data_handler.get_features(data_annotation_.gene_id)

    # features_ = Genomics.entries_for_gene_annotation(data_annotationn_, args.window, features_metadata)
    if weights:
        x_w = robjects.FloatVector(features_.w.values)
    else:
        x_w = None

    if features_.shape[0] == 0:
        logging.log(9, "No features available")
        return
    logging.log(8, "Loading genotype data")

    features_data_, features_ = features_handler.load_features(features_, [x for x in d_['individual']])
    # features_data_ = Parquet._read(features, [x for x in features_.id.values],
    #                                specific_individuals=[x for x in d_["individual"]])

#    assert features_.shape[0] == len(features_data_), "Features {}, data {}".format(features_.shape[0], len(features_data_))
#    loaded_features = [v for v in features_data_.keys()]
#    loaded_features.remove('individual')
#    features_ = features_.
    logging.log(8, "training")
    weights, summary = train(features_data_, features_, d_, data_annotation_, x_w, not args.dont_prune, nested_folds)

    if weights.shape[0] == 0:
        logging.log(9, "no weights, skipping")
        return

    logging.log(8, "saving")
    weights = weights.assign(gene=data_annotation_.gene_id). \
        merge(features_.rename(columns={"id": "feature", "allele_0": "ref_allele", "allele_1": "eff_allele"}), on="feature"). \
        rename(columns={"feature": "varID"}). \
        assign(gene=gene_id_)

    weights = weights[["gene", "rsid", "varID", "ref_allele", "eff_allele", "weight"]]
    if args.output_rsids:
        weights.loc[weights.rsid == "NA", "rsid"] = weights.loc[weights.rsid == "NA", "varID"]
    w.write(weights.to_csv(sep="\t", index=False, header=False, na_rep="NA").encode())

    summary = summary. \
        assign(gene=gene_id_, genename=data_annotation_.gene_name,
               gene_type=data_annotation_.gene_type). \
        rename(columns={"n_features": "n_snps_in_window", "n_features_in_model": "n.snps.in.model",
                        "zscore_pval": "pred.perf.pval", "rho_avg_squared": "pred.perf.R2",
                        "cv_converged":"nested_cv_converged"})
    summary["pred.perf.qval"] = None
    summary = summary[summary_fields]
    s.write(summary.to_csv(sep="\t", index=False, header=False, na_rep="NA").encode())

    var_ids = [x for x in weights.varID.values]
    cov = numpy.cov([features_data_[k] for k in var_ids], ddof=1)
    ids = [x for x in weights.rsid.values] if args.output_rsids else var_ids
    cov = matrices._flatten_matrix_data([(gene_id_, ids, cov)])
    for cov_ in cov:
        l = "{} {} {} {}\n".format(cov_[0], cov_[1], cov_[2], cov_[3]).encode()
        c.write(l)

########################################################################################################################


    # @staticmethod
    # def _get_dapg_weights(x_weights, id_whitelist=None, pre_parsed=False):
    #     if pre_parsed:
    #         weights_df_lst = []
    #         for i in x_weights:
    #             weights_df_lst.append(Miscellaneous.dapg_preparsed(i))
    #         df = pandas.concat(weights_df_lst)
    #         df.w = 1 - df.w
    #         return df
    #     else:
    #         if id_whitelist is None:
    #             raise ValueError("Either id_whitelist or pre_parsed must be specified")
    #         if x_weights[1] == "PIP":
    #             w = Miscellaneous.dapg_signals(x_weights[0], float(x_weights[2]), id_whitelist)
    #             w = w.rename(columns={"gene": "gene_id", "pip": "w", "variant_id": "id"})
    #             w.w = 1 - w.w  # Less penalty to the more probable snps
    #         else:
    #             raise RuntimeError("unsupported weights argument")
    #         return w



########################################################################################################################
def run(args):

    logging.info("Starting")
    Utilities.ensure_requisite_folders(args.output_prefix)

    wp = args.output_prefix + "_weights.txt.gz"
    sp = args.output_prefix + "_summary.txt.gz"
    cp = args.output_prefix + "_covariance.txt.gz"
    r = args.output_prefix + "_run.txt.gz"

    out_files = [wp, sp, cp, r]
    for i in out_files:
        Utilities.ensure_no_file(i)


    logging.info("Opening pheno data")
    d_handler = Parquet.PhenoDataHandler(args.data, sub_batches=args.sub_batches,
                                         sub_batch=args.sub_batch)
    # data = pq.ParquetFile(args.data)
    # available_data = {x for x in data.metadata.schema.names}

    if args.features_weights:
        logging.info("Loading weights")
        weights = get_weights(args.features_weights, None)
        d_handler.add_features_weights(weights)


    if args.preparsed_weights:
        logging.info("Loading preparsed weights")
        weights = get_dapg_preparsed(args.preparsed_weights)
        d_handler.add_features_weights(weights)
        # x_weights = get_weights(args.features_weights, pre_parsed=True)
        # whitelist = { v for v in x_weights.variant_id}

    f_handler = Parquet.MultiFileGenoHandler(args.features,
                                             args.features_annotation)
    logging.info("Loading geno metadata")
    features_metadata = f_handler.load_metadata()
    d_handler.add_features_metadata(features_metadata)


    # TODO: Re -integrate data_annotation
    # if args.data_annotation:
    #     logging.info("Loading data annotation")
    #     data_annotation = StudyUtilities.load_gene_annotation(args.data_annotation, args.chromosome, args.sub_batches, args.sub_batch)
    #     data_annotation = data_annotation[data_annotation.gene_id.isin(available_data)]
    #     if args.gene_whitelist:
    #         logging.info("Applying gene whitelist")
    #         data_annotation = data_annotation[data_annotation.gene_id.isin(set(args.gene_whitelist))]
    #     logging.info("Kept %i entries", data_annotation.shape[0])
    # else:
    #     data_annotation = None

    # logging.info("Opening features annotation")
    # if not args.chromosome:
    #     features_metadata = pq.read_table(args.features_annotation).to_pandas()
    # else:
    #     features_metadata = pq.ParquetFile(args.features_annotation).read_row_group(args.chromosome-1).to_pandas()

    # if args.chromosome and args.sub_batches and data_annotation:
    #     logging.info("Trimming variants")
    #     features_metadata = StudyUtilities.trim_variant_metadata_on_gene_annotation(features_metadata, data_annotation, args.window)


    # if args.rsid_whitelist:
    #     logging.info("Filtering features annotation")
    #     whitelist = TextFileTools.load_list(args.rsid_whitelist)
    #     whitelist = set(whitelist)
    #     features_metadata = features_metadata[features_metadata.rsid.isin(whitelist)]


    # else:
    #     x_weights = None

    # if data_annotation is None:
    #     d_ = pq.ParquetFile(args.data)
    #     gene_lst = d_.metadata.schema.names
    #     gene_lst.remove('individual')
    #     dd = {'gene_name': gene_lst, 'gene_id': gene_lst,
    #           'gene_type': ['NA'] * len(gene_lst)}
    #     data_annotation = pandas.DataFrame(dd)

#    logging.info("Opening features")
#    features = pq.ParquetFile(args.features)

    logging.info("Setting R seed")
    s = numpy.random.randint(1e8)
    set_seed(s)
    if args.run_tag:
        d = pandas.DataFrame({"run":[args.run_tag], "cv_seed":[s]})[["run", "cv_seed"]]
        Utilities.save_dataframe(d, r)

    WEIGHTS_FIELDS=["gene", "rsid", "varID", "ref_allele", "eff_allele", "weight"]
    SUMMARY_FIELDS=["gene", "genename", "gene_type", "alpha", "n_snps_in_window", "n.snps.in.model",
                    "test_R2_avg", "test_R2_sd", "cv_R2_avg", "cv_R2_sd", "in_sample_R2", "nested_cv_fisher_pval",
                    "nested_cv_converged", "rho_avg", "rho_se", "rho_zscore", "pred.perf.R2", "pred.perf.pval", "pred.perf.qval"]

    train = train_elastic_net_wrapper if args.mode == "elastic_net" else train_ols

    with gzip.open(wp, "w") as w:
        w.write(("\t".join(WEIGHTS_FIELDS) + "\n").encode())
        with gzip.open(sp, "w") as s:
            s.write(("\t".join(SUMMARY_FIELDS) + "\n").encode())
            with gzip.open(cp, "w") as c:
                c.write("GENE RSID1 RSID2 VALUE\n".encode())
                for i,data_annotation_ in enumerate(d_handler.data_annotation.itertuples()):
                    if args.MAX_M and  i>=args.MAX_M:
                        logging.info("Early abort")
                        break
                    logging.log(9, "processing %i/%i:%s", i+1, d_handler.data_annotation.shape[0], data_annotation_.gene_id)
                    if args.repeat:
                        for j in range(0, args.repeat):
                            logging.log(9, "%i-th reiteration", j)
                            process(w, s, c, d_handler, data_annotation_,
                                    f_handler, d_handler.send_weights,
                                    SUMMARY_FIELDS, train, j, args.nested_cv_folds)
                    else:
                        process(w, s, c, d_handler, data_annotation_,
                                f_handler, d_handler.send_weights,
                                SUMMARY_FIELDS, train,
                                nested_folds=args.nested_cv_folds)

    logging.info("Finished")


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser("Train Elastic Net prediction models from GLMNET")
    parser.add_argument("--features_weights", nargs="+")
    parser.add_argument("-features")
    parser.add_argument("-features_annotation")
    parser.add_argument("-data", help="Phenotype data, parquet format")
    parser.add_argument("-data_annotation")
    parser.add_argument("-window", type = int)
    parser.add_argument("--run_tag")
    parser.add_argument("--output_rsids", action="store_true")
    parser.add_argument("--chromosome", type = int)
    parser.add_argument("--sub_batches", type = int)
    parser.add_argument("--sub_batch", type =int)
    parser.add_argument("--rsid_whitelist")
    parser.add_argument("--MAX_M", type=int)
    parser.add_argument("--mode", default="elastic_net", help="'elastic_net' or 'ols'")
    parser.add_argument("--gene_whitelist", nargs="+", default=None)
    parser.add_argument("--preparsed_weights", help="Pre-parsed dapg weights")
    parser.add_argument("--dont_prune", action="store_true")
    parser.add_argument("-output_prefix")
    parser.add_argument("-parsimony", default=10, type=int)
    parser.add_argument("--repeat", default=None, type=int)
    parser.add_argument("--nested_cv_folds", default=5, type=int)
    args = parser.parse_args()
    Logging.configure_logging(args.parsimony)

    initialize()

    run(args)
