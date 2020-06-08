import logging
import os
from ...individual_data.Utilities import StudyBasedContext
from ...individual_data import Utilities as StudyUtilities
from ...miscellaneous import PandasHelpers
from ...file_formats import Parquet

class Context(StudyBasedContext):
    def __init__(self, fastenloc_exe, study, gene_annotation, window):
        super(Context, self).__init__(study, gene_annotation, window)
        self.fastenloc_exe = fastenloc_exe

    def get_available_genes(self):
        pass

def context_from_args(args):
    logging.info("Creating context")

