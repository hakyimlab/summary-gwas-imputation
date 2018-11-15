#!/usr/bin/env python
__author__ = "alvaro barbeira"

import re
import logging
import os
import shutil

from genomic_tools_lib import Logging, Utilities

def run(args):
    if not args.reentrant:
        if os.path.exists(args.output_folder):
            logging.info("Output path exists. Nope.")
            return

    Utilities.maybe_create_folder(args.output_folder)


    logging.info("Checking input folder")
    r = re.compile(args.rule)
    folders = [x for x in sorted(os.listdir(args.input_folder)) if r.search(x)]
    if args.exclude:
        folders = [x for x in folders if not x in {y for y in args.exclude}]
    names = {}
    for f in folders:
        name = r.search(f).group(1)
        if not name in names: names[name] = []
        names[name].append(os.path.join(args.input_folder, f))


    _f = shutil.move if args.move else shutil.copy
    for name in sorted(names):
        logging.info("Processing %s", name)
        output_folder = os.path.join(args.output_folder, name)
        Utilities.maybe_create_folder(output_folder)

        for input_folder in names[name]:
            logging.log(8, "Processing %s", input_folder)
            files = os.listdir(input_folder)
            for file in files:
                i = os.path.join(input_folder, file)
                o = os.path.join(output_folder, file)
                _f(i, o)
    logging.info("Finished collapse")

if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser("Convert model training format data to parquet format ")
    parser.add_argument("-input_folder", help="Folder where sub folders can be found")
    parser.add_argument("-rule", help="Regexp to group input folders")
    parser.add_argument("-output_folder", help="Destination folder where contets will be saved")
    parser.add_argument("--reentrant", help="Lenient, multiple-run mode", action="store_true")
    parser.add_argument("--exclude", help="Skip these folders", nargs="+")
    parser.add_argument("--move", help="Wether to move or copy files", action="store_true")
    parser.add_argument("-parsimony", help="Log parsimony level. 1 is everything being logged. 10 is only high level messages, above 10 will hardly log anything", default = "10")

    args = parser.parse_args()

    Logging.configure_logging(int(args.parsimony))

    run(args)