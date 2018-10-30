__author__ = "alvaro barbeira"

import  logging
from ...data_management import GTExMisc
from ...file_formats import Parquet
from .RunTorus import _Context

class Context(_Context):
    def __init__(self, input_eqtl,  input_gene_annotation, intermediate_folder, output_folder, torus_exe, delete_intermediate, variant_whitelist=None):
        self.input_eqtl = input_eqtl
        self.variant_whitelist = variant_whitelist
        self.input_gene_annotation = input_gene_annotation
        self.intermediate_folder = intermediate_folder
        self.output_folder = output_folder
        self.torus_exe = torus_exe
        self.delete_intermediate = delete_intermediate

    def get_input_eqtl(self): return  self.input_eqtl
    def get_variant_whitelist(self): return self.variant_whitelist
    def get_input_gene_annotation(self): return self.input_gene_annotation
    def get_intermediate_folder(self): return self.intermediate_folder
    def get_torus_exe(self): return self.torus_exe
    def get_output_folder(self): return self.output_folder
    def get_delete_intermediate(self): return self.delete_intermediate

def context_from_args(args):
    whitelist = None
    if args.snp_annotation_file:
        logging.info("Loading SNP annotation file for its list of variants")
        whitelist = GTExMisc.load_gtex_variants(args.snp_annotation_file, args.frequency_filter)
    elif args.snp_annotation_from_parquet_metadata:
        logging.info("Loading Parquet metadata for its list of variants")
        whitelist = Parquet.variants_from_metadata(args.snp_annotation_from_parquet_metadata, args.frequency_filter)
    else:
        raise RuntimeError("A metadata file is required at the time")

    context = Context(args.eqtl, args.gene_annotation_file,
                args.intermediate_folder, args.output_folder, args.torus_command,
                      (not args.keep_intermediate_folder), variant_whitelist=whitelist)
    return context