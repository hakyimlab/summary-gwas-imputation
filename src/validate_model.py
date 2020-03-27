
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
from timeit import default_timer as timer

import pyarrow.parquet as pq

from genomic_tools_lib import Utilities, Logging
from genomic_tools_lib.data_management import TextFileTools
from genomic_tools_lib.individual_data import Utilities as StudyUtilities
from genomic_tools_lib.file_formats import Parquet, Miscellaneous
from genomic_tools_lib.miscellaneous import matrices, Genomics, Math



def loss(X, y_true, weights):
    if X.shape[1] != weights.shape[0]:
        s = "X and weights have incompatible " \
            "dimensions: {} and {}".format(X.shape[1], weights.shape[0])
        raise ValueError(s)

    if X.shape[0] != y_true.shape[0]:
        s = "X and y have incompatible " \
            "lengths: {} and {}".format(X.shape[0], y_true.shape[0])
    y_pred = np.dot(X, weights)
    score = r2_score(y_true, y_pred)
    correlation = pearsonr(y_true, y_pred)[0]
    return (score, correlation)

def load_weights(fp):
    return pandas.read_csv(fp, sep = "\t")

def validate(pheno, f_handler, d_handler, individuals):
    """

    :param pheno: str. Pheno name
    :param f_handler: instance of Parquet.MultiFileGenoHandler
    :param d_handler: instance of Parquet.PhenoDataHandler
    """
    pheno_data = d_handler.load_pheno(pheno, to_pandas=True)
    pheno_data = pheno_data.set_index('individual')

    features_ = d_handler.get_features(pheno)

    geno_data = f_handler.load_features(features_, individuals, pandas=True)




def run(args):
    start = timer()

    out_pre = os.path.dirname(args.out_fp)
    Utilities.ensure_requisite_folders(out_pre)
    Utilities.ensure_no_file(args.out_fp)

    logging.info("Opening pheno data")
    d_handler = Parquet.PhenoDataHandler(args.data, args.sub_batches,
                                         args.sub_batch)
    individuals = d_handler.get_individuals()
    logging.log(9, "Phenotype data for {} individuals".format(len(individuals)))

    logging.info("Loading weights")
    weights = load_weights(args.model_weights)
    var_lst = list(weights.id)
    d_handler.add_features_weights(weights)



    logging.info("Opening geno metadata")
    f_handler = Parquet.MultiFileGenoHandler(args.features,
                                             args.features_metadata)
    features_metadata = f_handler.load_metadata(whitelist=var_lst)



    results_cols = ['pheno', 'r2', 'correlation']
    results_header = "\t".join(results_cols) + "\n"
    with open(args.out_fp, 'w') as f:
        f.write(results_header)

    n_phenos = len(d_handler.data_annotation)
    for i, pheno in enumerate(d_handler.data_annotation.gene_id.unique()):
        logging.log(9, "processing %i/%i:%s", i+1, n_phenos, pheno)

        results = validate(pheno, f_handler, d_handler, individuals)


        weights_df, pheno_series = data_generator.get_pheno_model(i, pheno_map.map_pheno(i))
        dm, rsids = data_generator.get_design_matrix(i)
        weights_df = weights_df.reindex(rsids)
        val_weights = np.array(weights_df['weights'])
        loss_r2, loss_corr = loss(dm, np.array(pheno_series), val_weights)
        with open(result_fp, 'a') as f:
            f.write("{}\t{}\t{}\n".format(i, loss_r2, loss_corr))
        # results_df = results_df.append({'r2': loss_val,
        #                                 'pheno': i})
        # logging.info("Writing validation results to {}".format(result_fp))
        # results_df.to_csv(result_fp, sep = "\t", index=False)
        logging.info("Wrote results to {}".format(result_fp))


if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser("Validate prediction models with held-out "
                                     "data.")

    parser.add_argument('-features', help="Genotype file")
    parser.add_argument('-features_annotation', help="Geno metadata file")
    parser.add_argument('-model_weights', help="Predicion model weights")
    parser.add_argument('-out_fp')
    parser.add_argument('-data', help="Phenotype data.")
    parser.add_argument('--sub_batches', type = int)
    parser.add_argument('--sub_batch', type = int)
    parser.add_argument('-parsimony', type = int, default = 8)

    args = parser.parse_args()

    Logging.configure_logging(level=args.parsimony, with_date=True)


    run(args)
