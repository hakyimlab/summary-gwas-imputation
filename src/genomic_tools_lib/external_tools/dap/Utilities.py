__author__ = "alvaro barbeira"

import logging
import os
from . import RunDAP
from ...individual_data.Utilities import StudyBasedContext
from ...individual_data import Utilities as StudyUtilities
from ...miscellaneous import PandasHelpers
from ...file_formats import Parquet

class Context(StudyBasedContext, RunDAP._Context):
    def __init__(self, dap_exe, grid_file, prior_folder, intermediate_folder,
                 output_folder, study, gene_annotation, window,
                 delete_intermediate, options, special_gene_handling):
        super().__init__(study, gene_annotation, window)
        self.dap_exe = dap_exe
        self.grid_file = grid_file
        self.prior_folder = prior_folder
        self.intermediate_folder = intermediate_folder
        self.output_folder = output_folder
        if self.prior_folder is not None:
            self.gene_to_prior = {os.path.splitext(x)[0]:os.path.join(prior_folder,x) for x in os.listdir(prior_folder)}
        else:
            self.gene_to_prior = None
        self.delete_intermediate = delete_intermediate
        self.options = options
        self.special_gene_handling = special_gene_handling
        if self.special_gene_handling:
            self.get_available_genes = self._available_genes_special_handling
        else:
            self.get_available_genes = self._available_genes


    def get_dap_exe(self): return self.dap_exe
    def get_grid_file_path(self): return self.grid_file
    def get_intermediate_folder(self): return self.intermediate_folder
    def get_output_folder(self): return self.output_folder
    def get_delete_intermediate(self):return self.delete_intermediate
    def get_options(self): return self.options
    def get_prior_file_path(self, gene):
        if self.gene_to_prior is not None:
            return self.gene_to_prior[gene]
        else:
            return None

    def _available_genes(self):
        if not self.available_genes_:
            g = StudyUtilities.available_genes(self.study, self.gene_annotation)
            if self.gene_to_prior is not None:
                self.available_genes_ = [x for x in g if x in self.gene_to_prior]
            else:
                self.available_genes_ = g
        return self.available_genes_

    def _available_genes_special_handling(self):
        if not self.available_genes_:
            gene_ids = self.gene_annotation.gene_id.values
            gene_names = self.gene_annotation.gene_name.values
            ll = len(gene_ids)
            self.available_genes_ = {gene_names[i]: gene_ids[i] for i in range(ll)}
        return self.available_genes_



def context_from_args(args):
    logging.info("Creating context")

    if args.gene_annotation is not None:
        gene_annotation = StudyUtilities.load_gene_annotation(args.gene_annotation, args.chromosome)
    else:
        gene_annotation = StudyUtilities.gene_annotation_from_parquet_phenotype(args.parquet_phenotype)
    if args.sub_batches and args.sub_batch is not None:
        logging.log(9, "Trimming gene annotation on sub-batches")
        gene_annotation = PandasHelpers.sub_batch(gene_annotation, args.sub_batches, args.sub_batch)

    if args.regions is not None and args.chromosome is not None:
        gene_annotation = StudyUtilities.expand_gene_annotation_from_regions(gene_annotation, args.regions, args.chromosome)

    logging.info("Loading study")
    p_ = (lambda x: StudyUtilities.trim_variant_metadata_on_gene_annotation(x, gene_annotation, args.window, _log_level_v=5)) if args.sub_batches else None
    study = Parquet.study_from_parquet(args.parquet_genotype, args.parquet_genotype_metadata, args.parquet_phenotype, args.parquet_covariate, post_process_variants_metadata=p_, frequency_filter=args.frequency_filter, chromosome=args.chromosome)
    special_gene_handling = True

    options= {}
    if args.options:
        options = {"-"+x[0]:x[1] for x in args.options}

    context = Context(args.dap_command, args.grid_file, args.priors_folder,
                      args.intermediate_folder, args.output_folder, study,
                      gene_annotation, args.window,
                      (not args.keep_intermediate_folder), options,
                      special_gene_handling)
    return context
