#!/usr/bin/env python
__author__ = "alvaro barbeira"

import os
import logging
from timeit import default_timer as timer

from genomic_tools_lib import Logging, Utilities
from genomic_tools_lib.file_formats import DBSnp

def run(args):
    if (not args.output and not args.output_blacklist) or (args.output and args.output_blacklist):
        logging.info("Provide only one output argument")
        return


    if args.output and os.path.exists(args.output):
        logging.info("Output path %s exists. Nope.", args.output)
        return

    if args.output_blacklist and os.path.exists(args.output_blacklist):
        logging.info("Output path for skipped variants %s exists. Nope.", args.output_blacklist)
        return

    start = timer()
    logging.info("Started parsing DB SNP file")

    if args.output:
        Utilities.ensure_requisite_folders(args.output)
        entries = Utilities.lineify(DBSnp.generate(args.input, args.fields, args.keep_zero_based, recode_observed=args.recode_observed))
        Utilities.write_iterable_to_file(entries, args.output, Utilities.to_line(args.fields))
    else:
        Utilities.ensure_requisite_folders(args.output_blacklist)
        entries = Utilities.lineify(DBSnp.generate_skips(args.input, args.fields, args.keep_zero_based, recode_observed=args.recode_observed))
        Utilities.write_iterable_to_file(entries, args.output_blacklist, Utilities.to_line(args.fields+ ["reason"]))

    end = timer()
    logging.info("Finished parsing at %s", str(end-start))

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser("parse/convert dbsnp file")
    parser.add_argument("-keep_zero_based", action="store_true")
    parser.add_argument("-input", help="input DB SNP file")
    parser.add_argument("-output", help="output DB SNP file")
    parser.add_argument("-output_blacklist", help="output variants to skip from DB SNP file")
    parser.add_argument("-fields", help="fields to extract", default=["chromosome", "start", "end", "name"], nargs="+")
    parser.add_argument("--recode_observed", action="store_true")
    parser.add_argument("-parsimony", help="log output parsimony", type=int, default=logging.INFO)
    args = parser.parse_args()

    Logging.configure_logging(args.parsimony)
    run(args)