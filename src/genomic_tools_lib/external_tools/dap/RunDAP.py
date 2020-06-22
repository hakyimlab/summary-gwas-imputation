__author__ = "alvaro barbeira"

import os
import re
import numpy
import shutil
import logging
import traceback
from collections import namedtuple
from subprocess import call, PIPE

from ...file_formats import SBAM
from ...individual_data.Utilities import _StudyBasedContext
from ...individual_data import Utilities as StudyUtilities
from ...Exceptions import ReportableException
from ... import Utilities

COMMAND_WITH_PRIORS = """#!/usr/bin/env bash

unset OMP_NUM_THREADS

[ -d {OUTPUT_DIR} ] || mkdir -p {OUTPUT_DIR}
[ -d {INTERMEDIATE_DIR} ] || mkdir -p {INTERMEDIATE_DIR}

{dap} \\
-d {data} \\
-prior {prior} \\
{extra} \\
> {output}
"""
COMMAND_NO_PRIORS = """#!/usr/bin/env bash

unset OMP_NUM_THREADS

[ -d {OUTPUT_DIR} ] || mkdir -p {OUTPUT_DIR}
[ -d {INTERMEDIATE_DIR} ] || mkdir -p {INTERMEDIATE_DIR}

{dap} \\
-d {data} \\
{extra} \\
-o {output}
"""
Stats = namedtuple("Stats", ["gene", "status"])

SEGFALT_RE = re.compile("[Ss]egmentation [Ff]ault|[Ss]egfault")

# TODO: abstract and refactor into utilities
class _Context(_StudyBasedContext):
    def get_dap_exe(self): raise RuntimeError("Not implemented")
    def get_grid_file_path(self): raise RuntimeError("Not implemented")
    def get_options(self): raise RuntimeError("Not implemented")
    def get_prior_file_path(self, gene): raise RuntimeError("Not implemented")
    def get_intermediate_folder(self): raise RuntimeError("Not implemented")
    def get_output_folder(self): raise RuntimeError("Not implemented")
    def get_delete_intermediate(self): raise RuntimeError("Not implemented")

def run_dap(context, gene_id, gene_name):
    stats = _stats(gene_name)
    try:
        _run_dap(context, gene_id, gene_name)
    except ReportableException as ex:
        status = Utilities.ERROR_REGEXP.sub('_', ex.msg)
        stats = _stats(gene_name, status=status)
        logging.info("Reportable exception running dap: %s", ex.msg)
    except Exception as ex:
        msg = '{0}'.format(type(ex))  # .replace('"', "'")
        status = '"{0}"'.format(msg)
        stats = _stats(gene_name,  status=status)
        logging.info("Exception running dap:\n%s", traceback.format_exc())
    finally:
        if context.get_delete_intermediate():
            folder = _intermediate_folder(context, gene_name)
            if os.path.exists(folder):
                shutil.rmtree(folder)

    return stats

def _run_dap(context, gene_id, gene_name):
    intermediate = _intermediate_folder(context, gene_name)
    os.makedirs(intermediate)
    study = StudyUtilities._get_study_for_gene(context, gene_id, gene_name, rename_pheno=None)

    SBAM.save_study(study, intermediate)

    run_dap_command(context, gene_id, gene_name)

r_ = re.compile(r"\\\n[\s]+\\\n")
def _render(s):
    while r_.search(s):
        s = r_.sub("\\\n", s) #substitute empty lines on missing values
    return s

def dap_command(context, gene_id, gene_name):
    args = {"dap":context.get_dap_exe(),
            "data": _study_path(context, gene_id, gene_name),
            "grid": context.get_grid_file_path(),
             "prior": context.get_prior_file_path(gene_name),
             "output": _output(context, gene_name),
             "OUTPUT_DIR": context.get_output_folder(),
             "INTERMEDIATE_DIR": context.get_intermediate_folder()}

    options = context.get_options()
    if len(options):
        args["extra"] = " ".join(["{} {}".format(k, str(v)) for k,v in options.items()])
    else:
        args['extra'] = ""

    if args['prior'] is not None:
        command = COMMAND_WITH_PRIORS.format(**args)
    else:
        command = COMMAND_NO_PRIORS.format(**args)
#The following is not currently supported in dao-g. Sup.
#-t 1 \\
#-it 0.05 > {output}
    command = _render(command)
    return command

def run_dap_command(context, gene_id, gene_name):
    command = dap_command(context, gene_id, gene_name)
    script_path = _script_path(context, gene_name)
    with open(script_path, "w") as script:
        script.write(command)

    _o = os.path.join(_intermediate_folder(context, gene_name), "dap.o")
    _e = os.path.join(_intermediate_folder(context, gene_name), "dap.e")
    with open(_o, "w") as o:
        with open(_e, "w") as e:
            result = call(["bash", script_path], stderr=e, stdout=o)

    with open(_e, "r") as err_file:
        if SEGFALT_RE.search(err_file.read()):
            raise ReportableException("Segmentation Fault")

def _stats(gene, status=numpy.nan):
    return Stats(gene=gene, status=status)

########################################################################################################################

def _intermediate_folder(context, gene): return os.path.join(context.get_intermediate_folder(), gene)
def _study_path(context, gene_id, gene_name): return os.path.join(_intermediate_folder(context, gene_name), gene_id+".txt")
def _script_path(context, gene): return os.path.join(_intermediate_folder(context, gene), gene+".sh")
def _output(context, gene): return os.path.join(context.get_output_folder(), gene+".dap.txt")

########################################################################################################################

def data_frame_from_stats(data):
    return Utilities.to_dataframe(data, list(Stats._fields))
