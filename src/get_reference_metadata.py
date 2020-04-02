#!/usr/bin/env python
__author__ = "alvaro barbeira"

import os
import logging
import pandas

from genomic_tools_lib import Utilities, Logging
from genomic_tools_lib.data_management import GTExMisc
from genomic_tools_lib.file_formats import ModelTraining
from genomic_tools_lib.individual_data import Genotype

def run(args):
    if os.path.exists(args.output):
        logging.info("Output exists. Nope")
        return

    filters = {x[0]:x[1:] for x in args.filter}

    maf_filter = float(filters["MAF"][0]) if "MAF" in filters else None
    logging.info("Loading GTEX variant map")
    gtex_snp_key = GTExMisc.load_gtex_variant_to_rsid(args.annotation[0], args.rsid_column)

    logging.info("Processing genotype")
    m = []
    for mean, metadata, ids in ModelTraining.dosage_generator(args.genotype, gtex_snp_key, dosage_conversion=ModelTraining._mean, do_none=True):
        if maf_filter:
            f = mean / 2 if mean < 1 else 1 - mean / 2
            if f<maf_filter:
                continue
        m.append(metadata)

    m = Utilities.to_dataframe(m, [x[1] for x in Genotype.MetadataTFE.order])
    if "TOP_CHR_POS_BY_FREQ" in filters:
        logging.info("Simplifying multi-allelic variants")
        m = Genotype._monoallelic_by_frequency(m)

    logging.info("Saving...")
    Utilities.save_dataframe(m, args.output)
    logging.info("Finished")


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser("Parse a genotype text file for its variant metadata")
    parser.add_argument("-genotype", help="Path to genotype file")
    parser.add_argument("-annotation", help="Annotation file", nargs="+")
    parser.add_argument("-output", help = "Where to save the file")
    parser.add_argument("-filter", help="What to apply", nargs="+", action="append")
    parser.add_argument("-parsimony", help="Log parsimony", type=int, default=10)
    parser.add_argument("-rsid_column", help = "Column name with rsid variant identifiers", default = "rs_id_dbSNP150_GRCh38p7")
    args = parser.parse_args()

    Logging.configure_logging(args.parsimony)

    run(args)

