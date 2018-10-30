#! /usr/bin/env python
__author__ = "alvaro barbeira"

import os

from genomic_tools_lib import Logging, Utilities
from genomic_tools_lib.individual_data import Simulate
from genomic_tools_lib.file_formats import BIMBAM, Parquet, SBAM

def save_study(study, selected_snps, simulated_gencode, prefix, _save):
    Utilities.ensure_requisite_folders(prefix)
    _save(study)

    selected_snps_ = prefix + ".selected_snps.txt.gz"
    Utilities.write_iterable_to_file(selected_snps, selected_snps_)

    gencode_path = os.path.join(os.path.split(prefix)[0],  "gene_annotation.txt.gz")
    Utilities.save_dataframe(simulated_gencode, gencode_path)

def run(args):
    if not (args.bimbam_output_prefix or  args.parquet_output_prefix or args.sbam_output_folder):
        raise RuntimeError("Need output argument")

    #reproducibility. Add argument for different seed.
    Simulate.reset_seed()

    study, selected_snps, gene_annotation = Simulate.simulate_bslmm_study(args.snps_per_chromosome)

    if args.bimbam_output_prefix:
        _save = lambda study: BIMBAM.save_study(study, args.bimbam_output_prefix)
        save_study(study, selected_snps, gene_annotation, args.bimbam_output_prefix, _save)
    if args.parquet_output_prefix:
        _save = lambda study: Parquet.save_study(study, args.parquet_output_prefix)
        save_study(study, selected_snps, gene_annotation, args.parquet_output_prefix, _save)
    if args.sbam_output_folder:
        _save = lambda study: SBAM.save_study(study, args.sbam_output_folder)
        save_study(study, selected_snps, gene_annotation, os.path.join(args.sbam_output_folder, "_"), _save)

if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser("Generate Simulated Study")
    parser.add_argument("-bimbam_output_prefix", help="prefix for BIMBAM output")
    parser.add_argument("-parquet_output_prefix", help="prefix for PARQUET output")
    parser.add_argument("-sbam_output_folder", help="prefix for PARQUET output")
    parser.add_argument("-snps_per_chromosome", help="How many snps simulate per chromosome", type=int)
    parser.add_argument("-verbosity", help="Log verbosity level. 1 is everything being logged. 10 is only high level messages, above 10 will hardly log anything", default = "10")

    args = parser.parse_args()

    Logging.configure_logging(int(args.verbosity))

    run(args)