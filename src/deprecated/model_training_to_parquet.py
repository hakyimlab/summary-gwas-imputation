#!/usr/bin/env python
__author__ = "alvaro barbeira"

import logging

from timeit import default_timer as timer

from genomic_tools_lib import Logging
from genomic_tools_lib import Utilities
from genomic_tools_lib.file_formats import ModelTraining, Parquet
from genomic_tools_lib.data_management import KeyedDataSource


def process_phenotype(path, name, output_prefix):
    pheno = ModelTraining.load_variable_file(path)
    pheno_path = output_prefix + ".expression." + name + ".parquet"
    Parquet.save_variable(pheno_path, pheno)

def run(args):
    start = timer()
    Utilities.ensure_requisite_folders(args.output_prefix)
    logging.info("Loading SNP annotation")
    snp_key = KeyedDataSource.load_data(args.snp_annotation_file, "varID", "rsid_dbSNP150", should_skip=KeyedDataSource.skip_na)

    logging.info("Loading Genotype")
    genotype, individual_ids = ModelTraining.load_genotype_folder(args.input_genotype_folder, args.input_genotype_file_pattern, snp_key)

    logging.info("Saving Genotype")
    path_variant = args.output_prefix + ".variants.parquet"
    Parquet.save_variants(path_variant, genotype, individual_ids)

    path_metadata_variant = args.output_prefix + ".variants_metadata.parquet"
    Parquet.save_metadata(path_metadata_variant, genotype)

    logging.info("Processing Expression Phenotype")
    expression_logic = Utilities.file_logic(args.input_phenotype_folder, args.input_phenotype_expression_pattern)
    for row in expression_logic.itertuples():
        logging.info("Phenotype: %s", row.name)
        process_phenotype(row.path, row.name, args.output_prefix)
    end = timer()
    logging.info("Finished in %s", str(end-start))

if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser("Convert model training format data to parquet format ")
    parser.add_argument("-snp_annotation_file", help="file with snp annotation")

    parser.add_argument("-input_genotype_folder", help="Folder where genotype files are")
    parser.add_argument("-input_genotype_file_pattern", help="Regular expression to select files")
    parser.add_argument("-input_genotype_file", help="Alternative to the previous: just give the file")

    parser.add_argument("-input_phenotype_folder", help="Folder where genotype files are")
    parser.add_argument("-input_phenotype_expression_pattern", help="Regular expression to select files")
    parser.add_argument("-column_key", help="key to use as column names when transposing. Defaults to first column")
    parser.add_argument("-output_prefix", help="file to write transposed output")
    parser.add_argument("-columns_at_a_time", help="How many columns transpose at a time", type=int, default=50)
    parser.add_argument("-parsimony", help="Log parsimony level. 1 is everything being logged. 10 is only high level messages, above 10 will hardly log anything", default = "10")

    args = parser.parse_args()

    Logging.configure_logging(int(args.parsimony))

    run(args)