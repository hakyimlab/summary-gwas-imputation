__author__ = "alvaro barbeira"

import os
import re
import numpy
import shutil
import logging
import traceback
from collections import namedtuple
from subprocess import call

from ...file_formats import SBAM
from ...individual_data.Utilities import _StudyBasedContext
from ...individual_data import Utilities as StudyUtilities
from ...Exceptions import ReportableException
from ... import Utilities

Stats = namedtuple("Stats", ["gene", "status"])

# TODO: abstract and refactor into utilities
class _Context(_StudyBasedContext):
    def get_dap_exe(self): raise RuntimeError("Not implemented")
    def get_grid_file_path(self): raise RuntimeError("Not implemented")
    def get_options(self): raise RuntimeError("Not implemented")
    def get_prior_file_path(self, gene): raise RuntimeError("Not implemented")
    def get_intermediate_folder(self): raise RuntimeError("Not implemented")
    def get_output_folder(self): raise RuntimeError("Not implemented")
    def get_delete_intermediate(self): raise RuntimeError("Not implemented")

def run_dap(context, gene):
    stats = _stats(gene)
    try:
        _run_dap(context, gene)
    except ReportableException as ex:
        status = Utilities.ERROR_REGEXP.sub('_', ex.msg)
        stats = _stats(gene, status=status)
        logging.info("Reportable exception running dap: %s", ex.msg)
    except Exception as ex:
        msg = '{0}'.format(type(ex))  # .replace('"', "'")
        status = '"{0}"'.format(msg)
        stats = _stats(gene,  status=status)
        logging.info("Exception running dap:\n%s", traceback.format_exc())
    finally:
        if context.get_delete_intermediate():
            folder = _intermediate_folder(context, gene)
            if os.path.exists(folder):
                shutil.rmtree(folder)

    return stats

def _run_dap(context, gene):
    intermediate = _intermediate_folder(context, gene)
    os.makedirs(intermediate)
    study = StudyUtilities._get_study_for_gene(context, gene, rename_pheno=None)

    SBAM.save_study(study, intermediate)

    run_dap_command(context, gene)

r_ = re.compile(r"\\\n[\s]+\\\n")
def _render(s):
    while r_.search(s):
        s = r_.sub("\\\n", s) #substitute empty lines on missing values
    return s

def dap_command(context, gene):
    args = {"dap":context.get_dap_exe(),
         "data":_study_path(context, gene),
         "grid": context.get_grid_file_path(),
         "prior": context.get_prior_file_path(gene),
         "output": _output(context, gene),
         "OUTPUT_DIR": context.get_output_folder(),
         "INTERMEDIATE_DIR": context.get_intermediate_folder()}

    options = context.get_options()
    if len(options):
        args["extra"] = " ".join(["{} {}".format(k, str(v)) for k,v in options.items()])

    command = \
"""#!/usr/bin/env bash

[ -d {OUTPUT_DIR} ] || mkdir -p {OUTPUT_DIR}
[ -d {INTERMEDIATE_DIR} ] || mkdir -p {INTERMEDIATE_DIR}

{dap} \\
-d {data} \\
-prior {prior} \\
{extra} \\
-t 1 > {output}
""".format(**args)
#The following is not currently supported in dao-g. Sup.
#-t 1 \\
#-it 0.05 > {output}
    command = _render(command)
    return command

def run_dap_command(context, gene):
    command = dap_command(context, gene)
    script_path = _script_path(context, gene)
    with open(script_path, "w") as script:
        script.write(command)

    _o = os.path.join(_intermediate_folder(context, gene), "dap.o")
    _e = os.path.join(_intermediate_folder(context, gene), "dap.e")
    with open(_o, "w") as o:
        with open(_e, "w") as e:
            call(["bash", script_path], stderr=e, stdout=o)

def _stats(gene, status=numpy.nan):
    return Stats(gene=gene, status=status)

########################################################################################################################

def _intermediate_folder(context, gene): return os.path.join(context.get_intermediate_folder(), gene)
def _study_path(context, gene): return os.path.join(_intermediate_folder(context, gene), gene+".txt")
def _script_path(context, gene): return os.path.join(_intermediate_folder(context, gene), gene+".sh")
def _output(context, gene): return os.path.join(context.get_output_folder(), gene+".dap.txt")

########################################################################################################################

def data_frame_from_stats(data):
    return Utilities.to_dataframe(data, list(Stats._fields))
