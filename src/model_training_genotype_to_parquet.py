    #!/usr/bin/env python
__author__ = "alvaro barbeira"

import logging
import numpy
import pandas
from timeit import default_timer as timer

from genomic_tools_lib import Logging
from genomic_tools_lib import Utilities
from genomic_tools_lib.file_formats import ModelTraining, Parquet
from genomic_tools_lib.data_management import KeyedDataSource
from genomic_tools_lib.individual_data import Genotype, Utilities as GenotypeUtilities


def _filter(dosage, metadata, individual_ids, biallelic_only, maf_threshold):
    if biallelic_only and GenotypeUtilities._biallelic_filter(dosage, metadata, individual_ids):
        return True

    if maf_threshold and GenotypeUtilities._maf_filter_min_threshold(dosage, metadata, individual_ids, maf_threshold):
        return True

    return False

def get_filter(args, variant_key):
    if args.only_in_key:
        logging.info("Filtering on annotation")
        return (lambda dosage, metadata, individual_ids: metadata[Genotype.MetadataTF.ID] not in variant_key)
    elif args.biallelic_only or args.filter_maf:
        logging.info("Filter based on arguments")
        return (lambda dosage, metadata, individual_ids: _filter(dosage, metadata, individual_ids, args.biallelic_only, args.filter_maf))
    return None

def generate_single_backend(args, variant_key):
    logging.info("Loading Genotype")
    dosage_conversion = GenotypeUtilities.impute_to_mean_conversion if args.impute_to_mean else None
    dosage_filter = get_filter(args, variant_key)
    genotype, individual_ids = ModelTraining.load_genotype_file(args.input_genotype_file, variant_key, dosage_conversion, dosage_filter)

    if args.simplify_individual_id:
        logging.info("simplifying individual id")
        individual_ids = [x.split("_")[0] for x in individual_ids]

    logging.info("Saving Genotype")
    path_variant = args.output_prefix + ".variants.parquet"
    Parquet.save_variants(path_variant, genotype, individual_ids)

    logging.info("Saving metadata")
    path_metadata_variant = args.output_prefix + ".variants_metadata.parquet"
    Parquet._save_metadata(path_metadata_variant, genotype.get_variants_metadata())

def generate_multi_backend(args, variant_key):
    logging.info("Processing Genotype")
    dosage_conversion = GenotypeUtilities.impute_to_mean_conversion if args.impute_to_mean else None
    dosage_filter = get_filter(args, variant_key)
    metadata=[]

    for genotype, individual_ids in ModelTraining.load_genotype_file_by_chromosome(args.input_genotype_file, variant_key, dosage_conversion, dosage_filter):
        if args.simplify_individual_id:
            logging.info("simplifying individual id")
            individual_ids = [x.split("_")[0] for x in individual_ids]

        _m = genotype.get_variants_metadata()

        metadata.append(_m)
        _chr=_m.chromosome.values[0]
        logging.log(9, "Processing {}".format(_chr))
        _o = args.output_prefix + ".chr{}".format(_chr) + ".variants.parquet"
        Parquet.save_variants(_o, genotype, individual_ids)

    logging.info("Saving metadata")
    metadata = pandas.concat(metadata)
    path_metadata_variant = args.output_prefix + ".variants_metadata.parquet"
    Parquet._save_metadata(path_metadata_variant, metadata)

def get_variant_key(args):
    v_ = lambda x: numpy.nan if x=="." else x
    snp = args.snp_annotation_file
    variant_key = None
    if len(snp) == 1:
        variant_key = KeyedDataSource.load_data(snp[0], "variant_id", args.rsid_column,
                                            value_conversion=v_, key_filter=GenotypeUtilities.is_biallelic_variant)
    elif len(snp) == 2:
        if snp[1] == "METADATA":
            variant_key = KeyedDataSource.load_data(snp[0], "id", "rsid", value_conversion=v_)

    if not variant_key:
        raise  RuntimeError("Need right info to process snp metadata")

    return variant_key


def run(args):
    start = timer()
    Utilities.ensure_requisite_folders(args.output_prefix)

    logging.info("Loading SNP annotation")
    #TODO: make more generic
    variant_key = get_variant_key(args)
    if args.split_by_chromosome:
        generate_multi_backend(args, variant_key)
    else:
        generate_single_backend(args, variant_key)

    end = timer()
    logging.info("Finished in %s", str(end - start))


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser("Convert model training format data to parquet format ")
    parser.add_argument("-snp_annotation_file", help="file with snp annotation", nargs="+")
    parser.add_argument("-input_genotype_file", help="Alternative to the previous: just give the file")
    parser.add_argument("--only_in_key", help="Keep only variants in the snp annotation", action="store_true", default=False)
    parser.add_argument("--biallelic_only", help="Keep only biallelic snps", action="store_true", default=False)
    parser.add_argument("--filter_maf", help="If provided, filter variants outside a certain MAF threshold", type=float)
    parser.add_argument("--impute_to_mean", help="Impute missing values to the mean", action="store_true", default=False)
    parser.add_argument("-output_prefix", help="Prefix to generate genotype output")
    parser.add_argument("--split_by_chromosome", help="Split the genotype in one file per chromosome", action="store_true")
    parser.add_argument("--simplify_individual_id", help="Split the genotype in one file per chromosome", action="store_true")
    parser.add_argument("-parsimony", help="Log parsimony level. 1 is everything being logged. 10 is only high level messages, above 10 will hardly log anything", default="10")
    parser.add_argument("-rsid_column", help = "Column name with rsid variant identifiers", default = "rs_id_dbSNP150_GRCh38p7")

    args = parser.parse_args()

    Logging.configure_logging(int(args.parsimony))

    run(args)
