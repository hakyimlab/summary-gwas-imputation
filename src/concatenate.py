__author__ = "alvaro barbeira"
from genomic_tools_lib import Utilities, Logging
import os
import re
import gzip
import logging

def _key(x, r, sort_groups):
    k = ()
    s = r.search(x)
    for g in sort_groups:
        k = k + ( int(s.group(g)),)
    return k

def run(args):
    r = re.compile(args.pattern)
    files = [x for x in os.listdir(args.folder) if r.search(x)]
    if args.sort_groups:
        files = sorted(files, key=lambda x: _key(x, r, args.sort_groups))

    output_firstline = True
    Utilities.ensure_requisite_folders(args.output)

    logging.info("Starting concatenation")
    with gzip.open(args.output, "w") as o:
        for file in files:
            path = os.path.join(args.folder, file)
            logging.log(9, "Opening %s", path)
            for i, line in Utilities.iterate_file(path):
                if i==0:
                    if output_firstline:
                        o.write(line.encode())
                        if not args.headerless:
                            output_firstline = False
                    continue
                o.write(line.encode())

    logging.info("Finished")


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser("Post summary imputation results")
    parser.add_argument("-folder", help="How far to extend in each direction when searching for variants")
    parser.add_argument("-pattern", help="Work only with one chromosome")
    parser.add_argument("-output", help="Where to save stuff")
    parser.add_argument("--headerless", help="concatenate all lines", action="store_true")
    parser.add_argument("--sort_groups", help="what to extract and sort for", nargs="+", type=int)
    parser.add_argument("-parsimony",
                        help="Log verbosity level. 1 is everything being logged. 10 is only high level messages, above 10 will hardly log anything",
                        default="10")

    args = parser.parse_args()

    Logging.configure_logging(int(args.parsimony))

    run(args)