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

def get_weights(x_weights, id_whitelist, pre_parsed=False):
    if pre_parsed:
        return Miscellaneous.dapg_preparsed(x_weights)
    else:
        if x_weights[1] == "PIP":
            w = Miscellaneous.dapg_signals(x_weights[0], float(x_weights[2]), id_whitelist)
            w = w.rename(columns={"gene":"gene_id", "pip":"w", "variant_id":"id"})
            w.w = 1 - w.w #Less penalty to the more probable snps
        else:
            raise RuntimeError("unsupported weights argument")
        return w

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

def process(w, s, c, data, data_annotation_, features_handler, features_metadata,
            features, x_weights, summary_fields, train, postfix=None,
            nested_folds=10):
    """

    :param w: weights file handle
    :param s: summary file handle
    :param c: covariance file handle
    :param data: ParquetFile phenotype data
    :param data_annotation_: Pandas tuple with attributes 'gene_name', 'gene_id', 'gene_type'
    :param features_handler: class FeaturesHandler
    :param features_metadata: Pandas DataFrame
    ???  :param features: ParquetFile
    :param x_weights:
    :param summary_fields: list of strings
    :param train: function. Training function
    :param postfix: int. If doing repeats, this is the repeat number.
    :param nested_folds: int. Number of nested folds.
    :return:
    """
    gene_id_ = data_annotation_.gene_id if postfix is None else "{}-{}".format(data_annotation_.gene_id, postfix)
    logging.log(8, "loading phenotype data")
    d_ = Parquet._read(data, [data_annotation_.gene_id])

    features_ = Genomics.entries_for_gene_annotation(data_annotation_, args.window, features_metadata)
    print(type(features))
    print(features.head())
    if x_weights is not None:
        x_w = features_[["id"]].merge(x_weights[x_weights.gene_id == data_annotation_.gene_id], on="id")
        features_ = features_[features_.id.isin(x_w.id)]
        x_w = robjects.FloatVector(x_w.w.values)
    else:
        x_w = None

    if features_.shape[0] == 0:
        logging.log(9, "No features available")
        return

    features_data_ = Parquet._read(features, [x for x in features_.id.values],
                                   specific_individuals=[x for x in d_["individual"]])

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
class FeaturesHandler:
    """
    This class is for loading parquet metadata and genotype files. Most of its
    functionality is meant to assist in the case that both metadata and genotype
    are split into 22 files, but it should be robust to the case where there is
    only one genotype or metadata file.
    """
    def __init__(self, features, metadata):
        """
        If either argument is a pattern for multiple files, it must be
        formattable with the argument 'chr'
        """
        self.m_features = self.check_if_formattable(features)
        if self.m_features:
            self.features = self.format_chrom_file_names(features)
        else:
            self.features = [features]

        self.m_metadata = self.check_if_formattable(metadata)
        if self.m_metadata:
            self.metadata = self.format_chrom_file_names(metadata)
        else:
            self.metadata = [metadata]

    @staticmethod
    def format_chrom_file_names(s):
        l = [s.format(chr=i) for i in range(1, 23)]
        return l
    @staticmethod
    def check_if_formattable(s):
        matches = re.findall('{(.*?)}', s)
        if len(matches) > 0 and matches[0] == 'chr':
            return True
        else:
            return False

    def load_metadata(self, whitelist=None):
        if whitelist is not None:
            df_lst = []
            for i in self.metadata:
                df_i = pq.read_table(i).to_pandas()
                if whitelist is not None:
                    df_i = df_i[df_i.rsid.isin(whitelist)]
                df_lst.append(df_i)
            return pandas.concat(df_lst)

    def load_features(self, metadata):
        """
        :param metadata: pandas DataFrame with columns 'variant_id' and
                'chromosome'
        :return:
        """
        if self.m_features:
            return self._load_features_multiple(metadata)
        else:
            return self._load_features_single(metadata)

    def _load_features_single(self, metadata):
        return pq.read_table(self.features[0],
                             columns=list(metadata.variant_id))

    def _load_features_multiple(self, metadata):
        for chr, group in metadata.groupby('chromosome'):
            chr_fp = self.features[chr - 1]
            chr_vars =  list(group.variant_id)
            chr_features = pq.read_table(chr_fp, columns=chr_vars)

        return

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


    logging.info("Opening data")
    data = pq.ParquetFile(args.data)
    available_data = {x for x in data.metadata.schema.names}

    if args.features_weights:
        logging.info("Loading weights")
        x_weights = get_weights(args.features_weights, pre_parsed=True)
        whitelist = { v for v in x_weights.variant_id}
    else:
        x_weights = None
        whitelist = None

    features_handler = FeaturesHandler(args.featues, args.features_metadata)
    features_metadata = features_handler.load_metadata(whitelist=whitelist)

    if args.data_annotation:
        logging.info("Loading data annotation")
        data_annotation = StudyUtilities.load_gene_annotation(args.data_annotation, args.chromosome, args.sub_batches, args.sub_batch)
        data_annotation = data_annotation[data_annotation.gene_id.isin(available_data)]
        if args.gene_whitelist:
            logging.info("Applying gene whitelist")
            data_annotation = data_annotation[data_annotation.gene_id.isin(set(args.gene_whitelist))]
        logging.info("Kept %i entries", data_annotation.shape[0])
    else:
        data_annotation = None

    # logging.info("Opening features annotation")
    # if not args.chromosome:
    #     features_metadata = pq.read_table(args.features_annotation).to_pandas()
    # else:
    #     features_metadata = pq.ParquetFile(args.features_annotation).read_row_group(args.chromosome-1).to_pandas()

    if args.chromosome and args.sub_batches and data_annotation:
        logging.info("Trimming variants")
        features_metadata = StudyUtilities.trim_variant_metadata_on_gene_annotation(features_metadata, data_annotation, args.window)

    # if args.rsid_whitelist:
    #     logging.info("Filtering features annotation")
    #     whitelist = TextFileTools.load_list(args.rsid_whitelist)
    #     whitelist = set(whitelist)
    #     features_metadata = features_metadata[features_metadata.rsid.isin(whitelist)]


    # else:
    #     x_weights = None

    if data_annotation is None:
        d_ = pq.ParquetFile(args.data)
        gene_lst = d_.metadata.schema.names
        gene_lst.remove('individual')
        dd = {'gene_name': gene_lst, 'gene_id': gene_lst,
              'gene_type': ['NA'] * len(gene_lst)}
        data_annotation = pandas.DataFrame(dd)

    logging.info("Opening features")
    features = pq.ParquetFile(args.features)

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
                for i,data_annotation_ in enumerate(data_annotation.itertuples()):
                    if args.MAX_M and  i>=args.MAX_M:
                        logging.info("Early abort")
                        break
                    logging.log(9, "processing %i/%i:%s", i+1, data_annotation.shape[0], data_annotation_.gene_id)
                    if args.repeat:
                        for j in range(0, args.repeat):
                            logging.log(9, "%i-th reiteration", j)
                            process(w, s, c, data, data_annotation_,
                                    features_handler, features_metadata,
                                    x_weights, SUMMARY_FIELDS, train, j,
                                    args.nested_cv_folds)
                    else:
                        process(w, s, c, data, data_annotation_,
                                features_handler, features_metadata, x_weights,
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
    parser.add_argument("--dont_prune", action="store_true")
    parser.add_argument("-output_prefix")
    parser.add_argument("-parsimony", default=10, type=int)
    parser.add_argument("--repeat", default=None, type=int)
    parser.add_argument("--nested_cv_folds", default=5, type=int)
    args = parser.parse_args()
    Logging.configure_logging(args.parsimony)

    initialize()

    run(args)
