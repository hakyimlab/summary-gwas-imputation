__author__ = "alvaro barbeira"

import sys
import os
import logging
import gzip
import numpy
import pandas
import collections

import pyarrow.parquet as pq

from genomic_tools_lib import Utilities, Logging
from genomic_tools_lib.data_management import TextFileTools
from genomic_tools_lib.individual_data import Utilities as StudyUtilities
from genomic_tools_lib.file_formats import Parquet, Miscellaneous
from genomic_tools_lib.miscellaneous import matrices, Genomics

########################################################################################################################
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
########################################################################################################################

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

def run(args):
    wp = args.output_prefix + "_weights.txt.gz"
    if os.path.exists(wp):
        logging.info("Weights output exists already, delete it or move it")
        return

    sp = args.output_prefix + "_summary.txt.gz"
    if os.path.exists(sp):
        logging.info("Summary output exists already, delete it or move it")
        return

    cp = args.output_prefix + "_covariance.txt.gz"
    if os.path.exists(wp):
        logging.info("covariance output exists already, delete it or move it")
        return

    r = args.output_prefix + "_run.txt.gz"
    if os.path.exists(wp):
        logging.info("run output exists already, delete it or move it")
        return

    logging.info("Starting")
    Utilities.ensure_requisite_folders(args.output_prefix)

    logging.info("Opening data")
    data = pq.ParquetFile(args.data)
    available_data = {x for x in data.metadata.schema.names}

    logging.info("Loading data annotation")
    data_annotation = StudyUtilities.load_gene_annotation(args.data_annotation, args.chromosome, args.sub_batches, args.sub_batch)
    data_annotation = data_annotation[data_annotation.gene_id.isin(available_data)]
    logging.info("Kept %i entries", data_annotation.shape[0])

    logging.info("Opening features annotation")
    if not args.chromosome:
        features_metadata = pq.read_table(args.features_annotation).to_pandas()
    else:
        features_metadata = pq.ParquetFile(args.features_annotation).read_row_group(args.chromosome-1).to_pandas()

    if args.chromosome and args.sub_batches:
        logging.info("Trimming variants")
        features_metadata = StudyUtilities.trim_variant_metadata_on_gene_annotation(features_metadata, data_annotation, args.window)

    if args.rsid_whitelist:
        logging.info("Filtering features annotation")
        whitelist = TextFileTools.load_list(args.rsid_whitelist)
        whitelist = set(whitelist)
        features_metadata = features_metadata[features_metadata.rsid.isin(whitelist)]

    if args.features_weights:
        logging.info("Loading weights")
        x_weights = get_weights(args.features_weights, {x for x in features_metadata.id})
        logging.info("Filtering features metadata to those available in weights")
        features_metadata = features_metadata[features_metadata.id.isin(x_weights.id)]
        logging.info("Kept %d entries", features_metadata.shape[0])
    else:
        x_weights = None

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
                          "rho_avg", "rho_se", "rho_zscore", "pred.perf.R2", "pred.perf.pval", "pred.perf.qval"]


    with gzip.open(wp, "w") as w:
        w.write(("\t".join(WEIGHTS_FIELDS) + "\n").encode())
        with gzip.open(sp, "w") as s:
            s.write(("\t".join(SUMMARY_FIELDS) + "\n").encode())
            with gzip.open(cp, "w") as c:
                c.write("GENE RSID1 RSID2 VALUE\n".encode())
                for i,data_annotation_ in enumerate(data_annotation.itertuples()):
                    logging.log(9, "processing %i:%s", i, data_annotation_.gene_id)
                    logging.log(8, "loading data")
                    d_ = Parquet._read(data, [data_annotation_.gene_id])
                    features_ = Genomics.entries_for_gene_annotation(data_annotation_, args.window, features_metadata)

                    if x_weights is not None:
                        x_w = features_[["id"]].merge(x_weights[x_weights.gene_id == data_annotation_.gene_id], on="id")
                        features_ = features_[features_.id.isin(x_w.id)]
                        x_w = robjects.FloatVector(x_w.w.values)
                    else:
                        x_w = None

                    if features_.shape[0] == 0:
                        logging.log(9, "No features available")
                        continue

                    features_data_ = Parquet._read(features, [x for x in features_.id.values],
                                                   specific_individuals=[x for x in d_["individual"]])

                    logging.log(8, "training")
                    x = numpy.array([features_data_[v] for v in features_.id.values])
                    dimnames = robjects.ListVector([(1,robjects.StrVector(d_["individual"])), (2,robjects.StrVector(features_.id.values))])
                    x = robjects.r["matrix"](robjects.FloatVector(x.flatten()), ncol=features_.shape[0], dimnames=dimnames)
                    y = robjects.FloatVector(d_[data_annotation_.gene_id])
                    res = train_elastic_net(y, x, penalty_factor=x_w) #observation weights, not explanatory variable weight :( , x_weight = x_w)

                    weights = pandas2ri.ri2py(res[0])
                    if weights.shape[0] == 0:
                        logging.log(9, "no weights, skipping")
                        continue

                    # _save(d_, features_, features_data_, data_annotation_.gene_id)
                    # from IPython import embed; embed(); quit()

                    logging.log(8, "saving")
                    weights = weights.assign(gene= data_annotation_.gene_id).\
                        merge(features_.rename(columns={"id":"feature", "allele_0":"ref_allele", "allele_1":"eff_allele"}), on="feature").\
                        rename(columns={"feature":"varID"})

                    weights = weights[["gene", "rsid", "varID", "ref_allele", "eff_allele", "weight"]]
                    if args.output_rsids:
                        weights.loc[weights.rsid == "NA", "rsid"] = weights.loc[weights.rsid == "NA", "varID"]
                    w.write(weights.to_csv(sep="\t", index=False, header=False, na_rep="NA").encode())

                    summary = pandas2ri.ri2py(res[1]).\
                        assign(gene=data_annotation_.gene_id, genename=data_annotation_.gene_name, gene_type=data_annotation_.gene_type).\
                        rename(columns={"n_features":"n_snps_in_window", "n_features_in_model":"n.snps.in.model", "zscore_pval":"pred.perf.pval", "rho_avg_squared":"pred.perf.R2"})
                    summary["pred.perf.qval"] = None
                    summary = summary[SUMMARY_FIELDS]
                    s.write(summary.to_csv(sep="\t", index=False, header=False, na_rep="NA").encode())

                    var_ids = [x for x in weights.varID.values]
                    cov = numpy.cov([features_data_[k] for k in var_ids], ddof=1)
                    ids = [x for x in weights.rsid.values] if args.output_rsids else var_ids
                    cov = matrices._flatten_matrix_data([(data_annotation_.gene_id, ids, cov)])
                    for cov_ in cov:
                        l = "{} {} {} {}\n".format(cov_[0], cov_[1], cov_[2], cov_[3]).encode()
                        c.write(l)

                    write_header=False
                    
                    if args.MAX_M and  i>=args.MAX_M:
                        logging.info("Early abort")
                        break

    logging.info("Finished")


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser("Train Elastic Net prediction models from GLMNET")
    parser.add_argument("--features_weights", nargs="+")
    parser.add_argument("-features")
    parser.add_argument("-features_annotation")
    parser.add_argument("-data")
    parser.add_argument("-data_annotation")
    parser.add_argument("-window", type = int)
    parser.add_argument("--run_tag")
    parser.add_argument("--output_rsids", action="store_true")
    parser.add_argument("--chromosome", type = int)
    parser.add_argument("--sub_batches", type = int)
    parser.add_argument("--sub_batch", type =int)
    parser.add_argument("--rsid_whitelist")
    parser.add_argument("--MAX_M", type=int)
    parser.add_argument("-output_prefix")
    parser.add_argument("-parsimony", default=10, type=int)
    args = parser.parse_args()
    Logging.configure_logging(args.parsimony)

    initialize()

    run(args)