#!/usr/bin/env python3
__author__ = "alvaro barbeira"
import logging
import gzip
import pyliftover

from genomic_tools_lib import Utilities, Logging

def _l(liftover, chr, pos):
    _new_chromosome = "NA"
    _new_position = "NA"
    try:
        pos = int(pos)
        l_ = liftover.convert_coordinate(chr, pos)
        if l_:
            if len(l_) > 1:
                logging.warning("Liftover with more than one candidate: %s", t.variant_id)
            _new_chromosome = l_[0][0]
            _new_position = int(l_[0][1])
    except:
        pass
    return _new_chromosome, _new_position

def run(args):
    Utilities.ensure_requisite_folders(args.output)

    logging.info("starting lifting over.")
    liftover = pyliftover.LiftOver(args.liftover)
    with gzip.open(args.output, "w") as _o:
        with open(args.input) as _i:
            for i,line in enumerate(_i):
                if i ==0:
                    line = "\t".join(line.strip().split()) + "\n"
                    _o.write(line.encode())
                    continue

                try:
                    comps = line.strip().split()
                    chr = comps[0]
                    start = int(comps[1])
                    end = int(comps[2])

                    _chrs, _s = _l(liftover, chr, start)
                    _chre, _e = _l(liftover, chr, end)
                    if _chrs != _chre:
                        logging.warning("{}:{}:{} have different target chromosomes: {}/{}".format(chr, start, end, _chrs, _chre))
                    line = "{}\n".format("\t".join([_chrs, str(_s), str(_e)]))
                    _o.write(line.encode())
                except Exception as e:
                    logging.info("Error for: %s", line)


    logging.info("Finished lifting over.")

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser("Liftover ld regions file")
    parser.add_argument("-liftover", help="File with liftover chain")
    parser.add_argument("-input", help="region file")
    parser.add_argument("-output", help="Where the output should go")
    parser.add_argument("-verbosity", help="Log verbosity level. 1 is everything being logged. 10 is only high level messages, above 10 will hardly log anything", default = "10")

    args = parser.parse_args()

    Logging.configure_logging(int(args.verbosity))

    run(args)