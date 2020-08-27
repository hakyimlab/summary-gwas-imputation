import gzip
import pandas
import numpy
import logging
import re
import os

from genomic_tools_lib import Utilities, Logging
from genomic_tools_lib.miscellaneous import Models

__author__ = "alvaro barbeira"

def _read(folder, pattern, check=False):
    files = sorted([x for x in os.listdir(folder) if pattern.search(x)])
    files = [os.path.join(folder, x) for x in files]
    r =[]
    for file in files:
        try:
            d = pandas.read_table(file)
            if check:
                d = d.drop_duplicates()
            r.append(d)
        except:
            logging.info("issue opening file %s", file)
    return pandas.concat(r)

def _read_2(input_prefix, stem, check=False):
    path_ = os.path.split(input_prefix)
    r = re.compile(path_[1] + stem)
    return _read(path_[0], r, check=check)

def _sample_info(args):
    sample_info = None
    if args.sample_info:
        sample_info = pandas.DataFrame({"n_samples": [int(args.sample_info[0])], "population": [args.sample_info[1]], "tissue": [args.sample_info[2]]})
    return sample_info

def run(args):
    logging.info("Loading model summaries")
    extra = _read_2(args.input_prefix, "_summary.txt.gz")
    logging.log(7, "Found {} model summaries".format(extra.shape[0]))
    extra = extra[extra["n.snps.in.model"] > 0]
    if not args.skip_check:
        if "rho_avg" in extra:
            extra = extra[(extra["pred.perf.pval"] < 0.05) & (extra.rho_avg >0.1)]
        else:
            extra = extra[(extra["pred.perf.pval"] < 0.05)]
            extra = extra.assign(rho_avg = None)
        if not "pred.perf.qval" in extra:
            extra["pred.perf.qval"] = None

        if "nested_cv_converged" in extra:
            extra.nested_cv_converged = extra.nested_cv_converged.astype(numpy.int32)
    logging.log(7, "After checks, {} summaries remaining".format(extra.shape[0]))

    logging.info("Loading weights")
    weights = _read_2(args.input_prefix, "_weights.txt.gz")
    weights = weights[weights.gene.isin(extra.gene)]
    if args.coalesce_validation:
        logging.info("Loading validation scores")
        validation = _read_2(args.input_prefix, "_validation.txt.gz", check=False)
        logging.log(9, "Validation results for {} genes".format(validation.shape[0]))
        validation = validation[validation.gene.isin(extra.gene)]
        logging.log(5, "Validation results for {} genes".format(validation.shape[0]))
        Utilities.save_dataframe(validation,
                                 args.output_prefix + "_validation.txt")

    sample_info = _sample_info(args)
    if args.output_prefix:
        logging.info("Saving dbs and covariance")
        db = args.output_prefix + ".db"
        logging.info("Saving db")
        Models.create_model_db(db, extra, weights, sample_info)

        logging.info("Processing covariances")
        genes = {x for x in extra.gene}

        path_ = os.path.split(args.input_prefix)
        r = re.compile(path_[1] + "_covariance.txt.gz")
        files = sorted([x for x in os.listdir(path_[0]) if r.search(x)])
        files = [os.path.join(path_[0], x) for x in files]
        cov = args.output_prefix + ".txt.gz"
        with gzip.open(cov, "w") as cov_:
            cov_.write("GENE RSID1 RSID2 VALUE\n".encode())
            for nf,f in enumerate(files):
                logging.log(9, "file %i/%i: %s", nf, len(files), f)
                with gzip.open(f) as f_:
                    f_.readline()
                    for l in f_:
                        gene = l.decode().strip().split()[0]
                        if not gene in genes:
                            continue
                        cov_.write(l)

    if args.output_prefix_text:
        logging.info("Saving text output")
        Utilities.save_dataframe(weights, args.output_prefix_text+ "_t_weights.txt")
        Utilities.save_dataframe(extra, args.output_prefix_text + "_t_extra.txt")
    logging.info("Done")

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser("Compile model database and covariance from model training output")
    parser.add_argument("-input_prefix")
    parser.add_argument("-output_prefix")
    parser.add_argument("-output_prefix_text")
    parser.add_argument("-coalesce_validation", default=False, action='store_true')
    parser.add_argument("--sample_info", nargs="+", default=[])
    parser.add_argument("-parsimony", type=int, default=logging.INFO)
    parser.add_argument("--skip_check", default=False, action='store_true')

    args = parser.parse_args()
    Logging.configure_logging(args.parsimony)
    run(args)
