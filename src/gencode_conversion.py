#!/usr/bin/env python
import os
import logging
from genomic_tools_lib import Utilities, Logging
from genomic_tools_lib.file_formats import Gencode

__author__ = "alvaro barbeira"

#Quick and dirty hack for format needed elsewhere.
#TODO: formalize
def _reformat(gencode):
    gencode = gencode.rename(columns={"chromosome": "chr", "start_location": "start", "end_location": "end"})
    gencode.chr = gencode.chr.str.split("chr").str.get(1)
    return gencode[["chr", "gene_id", "gene_name", "start", "end", "gene_type"]]

def run(args):
    if os.path.exists(args.output):
        logging.info("Output exists. Nope.")
        return

    if args.output_column_map:
        selected = [x[0] for x in args.output_column_map]
    else:
        selected = [Gencode.GFTF.K_GENE_ID, Gencode.GFTF.K_GENE_NAME, Gencode.GFTF.K_GENE_TYPE]

    logging.info("Loading Gencode")
    gencode = Gencode.load(args.gencode_file,
        feature_type_whitelist={x for x in args.feature_type_whitelist},
        gene_type_white_list={x for x in args.gene_type_whitelist},
        transcript_type_whitelist={x for x in args.transcript_type_whitelist},
        selected_key_value_pairs=selected)
    #gencode = _reformat(gencode)
    logging.info("Converting format")
    if args.output_column_map:
        gencode = gencode.rename(columns={x[0]:x[1] for x in args.output_column_map})
        if "gene_version" in gencode and "gene_id" in gencode:
            gencode["gene_id"] = gencode.gene_id+ "." + gencode.gene_version
            keep = ["chromosome", "start_location", "end_location", "feature_type", "strand"]+[x[1] for x in args.output_column_map if x[1] not in {"gene_version"}]
            gencode = gencode[keep]
        else:
            gencode = gencode[["chromosome", "start_location", "end_location", "feature_type", "strand"] + [x[1] for x in
                                                                                                  args.output_column_map]]
    logging.info("Saving")
    Utilities.save_dataframe(gencode, args.output)
    logging.info("Finished")

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser("Convert gencode file into a table file")
    parser.add_argument("-gencode_file", help="Where to load file from")
    parser.add_argument("-output", help="Where to save")
    parser.add_argument("-gene_type_whitelist", help="Which types of genes to keep", default=[], nargs="+")
    parser.add_argument("-feature_type_whitelist", help="Which types of genes to keep", default=[], nargs="+")
    parser.add_argument("-transcript_type_whitelist", help="Which types of transcripts to keep", default=[], nargs="+")
    parser.add_argument("-output_column_map", help="Specify multiple key-value pairs to specify format conversion", nargs=2, action="append", default=[])
    parser.add_argument("-verbosity", help="Logging verbosity (actually loquacity)", type=int, default=10)
    args = parser.parse_args()

    Logging.configure_logging(args.verbosity)
    run(args)
