__author__ = "alvaro barbeira"

import logging
import os

import shutil
import traceback
from collections import namedtuple
from subprocess import call

import numpy
import pandas

from ...Exceptions import ReportableException
from ... import Utilities
from ...file_formats import BIMBAM
from ...individual_data.Utilities import _StudyBasedContext
from ...miscellaneous import matrices
from ...individual_data import Utilities as StudyUtilities

########################################################################################################################

class ResultStatsTF:
    GENE=0
    GENE_NAME=1
    N_SNPS_IN_MODEL=2
    STATUS = 3

    K_GENE="gene"
    K_GENE_NAME="gene_name"
    K_N_SNPS_IN_MODEL="n_snps_in_model"
    K_STATUS="status"

    order = [(GENE,K_GENE), (GENE_NAME, K_GENE_NAME), (N_SNPS_IN_MODEL, K_N_SNPS_IN_MODEL), (STATUS, K_STATUS)]

HyperParameters = namedtuple("HyperParameters",
    ["gene",
     "h_0_500", "h_0_025", "h_0_975",
     "pve_0_500", "pve_0_025", "pve_0_975",
     "rho_0_500", "rho_0_025", "rho_0_975",
     "pge_0_500", "pge_0_025", "pge_0_975",
     "pi_0_500", "pi_0_025", "pi_0_975",
     "n_gamma_0_500", "n_gamma_0_025", "n_gamma_0_975"
])

class _Context(_StudyBasedContext):
    def get_gemma_path(self): raise RuntimeError("Not implemented")
    def get_intermediate_path(self): raise RuntimeError("Not implemented")

########################################################################################################################

def run_gemma(context, gene):
    annotation = context.get_gene_annotation(gene)
    gene_name = annotation.gene_name
    weights, covariance_data, hyperparameters, stats = pandas.DataFrame(), matrices.matrix_data(gene, [], numpy.matrix([])), _hyper(gene), _stats(gene, gene_name)
    #weights, stats, covariance_data = _run_gemma(context, gene, annotation)
    try:
        #pass
        weights, covariance_data, hyperparameters, stats = _run_gemma(context, gene, annotation)
    except ReportableException as ex:
        status = Utilities.ERROR_REGEXP.sub('_', ex.msg)
        stats = _stats(gene, annotation.gene_name, status=status)
        logging.info("Reportable exception running gemma: %s", ex.msg)
    except Exception as ex:
        msg = '{0}'.format(type(ex))#.replace('"', "'")
        status = '"{0}"'.format(msg)
        stats = _stats(gene, annotation.gene_name, status=status)
        logging.info("Exception running gemma:\n%s", traceback.format_exc())
    finally:
        folder = _intermediate_folder(context, gene)
        if os.path.exists(folder):
            shutil.rmtree(folder)

    return weights, covariance_data, hyperparameters, stats

########################################################################################################################

def _intermediate_folder(context, gene):
    return os.path.join(context.get_intermediate_path(), gene)

def _study_prefix(context, gene):
    _i = _intermediate_folder(context, gene)
    return os.path.join(_i, "study")

def _output_folder(context, gene):
    _i = _intermediate_folder(context, gene)
    return os.path.join(_i, "results")

def _result(context, gene):
    return "results"

def _output_prefix(context, gene):
    return os.path.join(_output_folder(context, gene), _result(context, gene))

def _output_params(context, gene):
    return _output_prefix(context, gene) + ".param.txt"

def _output_hyperparams(context, gene):
    return _output_prefix(context, gene) + ".hyp.txt"

########################################################################################################################

def _stats(gene, gene_name, n_snps_in_model=numpy.nan, status=None):
    return (gene, gene_name, n_snps_in_model, status,)

def _hyper(gene):
    t = tuple([gene] + [None]*(len(HyperParameters._fields)-1))
    return HyperParameters._make(t)

def _run_gemma(context, gene, annotation):
    os.makedirs(_intermediate_folder(context, gene))
    study = StudyUtilities._get_study_for_gene(context, gene)
    study_prefix = _study_prefix(context, gene)
    BIMBAM.save_study(study, study_prefix)

    _run_gemma_exe(context, gene)

    weights = _get_weights(context, gene, study)
    covariance = numpy.cov(study.genotype.variants)
    covariance_data = matrices.matrix_data(gene, weights.rsid, covariance)
    hyper = _get_hyper(context, gene)
    stats = _stats(gene, annotation.gene_name, weights.shape[0])
    return weights, covariance_data, hyper, stats

def _get_weights(context, gene, study):
    params_path = _output_params(context, gene)
    if not os.path.exists(params_path):
        raise ReportableException("missing_gemma_parameters")
    params = pandas.read_table(params_path)
    params["weight"] = params.alpha + params.beta*params.gamma
    params = params[["rs", "weight"]].rename(columns={"rs":"rsid"}).assign(gene=gene)
    m = study.get_variants_metadata(params.rsid)
    m = m[["id", "allele_0", "allele_1"]].rename(columns={"id": "rsid", "allele_0":"ref_allele", "allele_1":"eff_allele"})
    w = m.merge(params, on="rsid")[["rsid", "gene", "weight", "ref_allele", "eff_allele"]]
    return w

def _get_hyper(context, gene):
    hyp_path = _output_hyperparams(context, gene)
    if not os.path.exists(hyp_path):
        raise ReportableException("missing_gemma_hyperparameters")
    #Take the second half of the values, assuming they are closer to Montecarlo Chain convergence
    hyp = pandas.read_table(hyp_path, sep="\s+")
    f_ = int(numpy.ceil(hyp.shape[0]/2))
    hyp = hyp.iloc[f_:,]

    r_ = (gene,)
    columns=["h", "pve", "rho", "pge", "pi", "n_gamma"]
    for k in columns:
        r_ += tuple(numpy.percentile(hyp[k], [50.0, 2.5, 97.5]))
    r_ = HyperParameters._make(r_)
    return r_

def _run_gemma_exe(context, gene):
    _intermediate = _intermediate_folder(context, gene)
    _script = os.path.join(_intermediate, "run.sh")
    _command = gemma_command(context, gene)
    with open(_script, "w") as f: f.write(_command)
    logging.log(6,"command:\n%s\n", _command)
    _o = os.path.join(_intermediate, "gemma_o.txt")
    _e = os.path.join(_intermediate, "gemma_e.txt")
    with open(_o, "w") as o:
        with open(_e, "w") as e:
            logging.log(8, "Running gemma for %s", gene)
            call(["bash", _script], stdout=o, stderr=e)

def gemma_command(context, gene):
    output_folder = _output_folder(context, gene)
    results = _result(context, gene)
    study = _study_prefix(context, gene)
    gemma = context.get_gemma_path()
    command = \
"""
#!/usr/bin/env bash
GEMMA_PATH={}
STUDY_PREFIX="{}"
OUTPUT_DIR="{}"
OUTPUT="{}"

#GEMMA has "output" folder as default
[ -d "$OUTPUT_DIR" ] || mkdir -p "$OUTPUT_DIR"

$GEMMA_PATH \\
-g "$STUDY_PREFIX.geno.txt.gz" \\
-p "$STUDY_PREFIX.pheno.txt" \\
-a "$STUDY_PREFIX.snp.txt" \\
-bslmm 1 \\
-outdir $OUTPUT_DIR \\
-o "$OUTPUT"
""".format(gemma, study, output_folder, results)
    return command

def dataframe_from_stats(stats):
    r = Utilities.to_dataframe(stats, [x[1] for x in ResultStatsTF.order])
    return r

def dataframe_from_covariance_data(data):
    d = matrices.matrices_data_to_dataframe(data, ["GENE", "RSID1", "RSID2", "VALUE"])
    return d

def dataframe_from_hyperparameters(data):
    d = Utilities.to_dataframe(data, list(HyperParameters._fields))
    return d