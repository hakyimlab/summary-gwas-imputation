__author__ = "alvaro barbeira"
import os
import logging
import pandas
from genomic_tools_lib import Logging, Utilities

def run(args):
    if os.path.exists(args.output_file):
        logging.info("Output %s exists. Nope", args.output_file)
        return

    results_order = []
    results = {}
    logging.info("Streaming file for groups")
    for i,line in Utilities.iterate_file(args.input_file):
        if i==0: continue

        comps = line.strip().split()
        key = comps[0]
        if not key in results:
            results_order.append(key)
            results[key] = 0
            logging.log(9, "Key: %s", str(key))
        results[key] += 1

    r = []
    logging.info("Producing output")
    for key in results_order:
        r.append((key, results[key]))
    r = pandas.DataFrame(r, columns=["key","count"])

    logging.info("Saving")
    Utilities.ensure_requisite_folders(args.output_file)
    Utilities.save_dataframe(r, args.output_file)

    logging.info("Finished.")

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser("Stream a file and count entries contiguous with a same key")
    parser.add_argument("-input_file")
    parser.add_argument("-output_file", help="Parquet Genotype variant metadata file")
    parser.add_argument("-parsimony", help="Log verbosity level. 1 is everything being logged. 10 is only high level messages, above 10 will hardly log anything", default = "10")

    args = parser.parse_args()

    Logging.configure_logging(int(args.parsimony))

    run(args)