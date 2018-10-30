__author__ = "alvaro barbeira"

import logging

from . import  RunGEMMA
from ...file_formats import Parquet
from ...miscellaneous import PandasHelpers
from ...individual_data.Utilities import StudyBasedContext
from ...individual_data import Utilities as StudyUtilities

class Context(StudyBasedContext, RunGEMMA._Context):
    def __init__(self, study, gemma, intermediate, gene_annotation, window):
        super().__init__(study, gene_annotation, window)
        self.gemma = gemma
        self.intermediate = intermediate

    def get_gemma_path(self): return self.gemma
    def get_intermediate_path(self):  return self.intermediate

def context_from_args(args):
    logging.info("Creating context")

    gene_annotation = StudyUtilities.load_gene_annotation(args.gene_annotation)
    if args.sub_batches and args.sub_batch is not None:
        logging.log(9, "Trimming gene annotation on sub-batches")
        gene_annotation = PandasHelpers.sub_batch(gene_annotation, args.sub_batches, args.sub_batch)

    logging.info("Loading study")
    p_ = (lambda x: StudyUtilities.trim_variant_metadata_on_gene_annotation(x, gene_annotation, args.window)) if args.sub_batches else None
    study = Parquet.study_from_parquet(args.parquet_genotype, args.parquet_genotype_metadata, args.parquet_phenotype, args.parquet_covariate, post_process_variants_metadata=p_)

    context = Context(study, args.gemma_command, args.intermediate_folder, gene_annotation, args.window)
    return context