#!/usr/bin/env python
__author__ = "alvaro barbeira"

from genomic_tools_lib import Utilities

def _to_line(comps, indexes, delimiter):
    return delimiter.join([comps[i] for i in indexes])

def _get_lines(input_file, columns, in_delim, out_delim):
    with Utilities.open_any(input_file) as _input:
        header = _input.readline().strip().split(in_delim)
        indexes = [header.index(x) for x in columns]
        yield _to_line(header, indexes, out_delim)
        for i, line in Utilities._iterate_file(_input):
            comps = line.strip().split()
            yield _to_line(comps, indexes, out_delim)

def run(args):
    delim = args.delimiter
    out_delim = delim if delim else "\t"
    Utilities.write_iterable_to_file(_get_lines(args.input_file, args.columns, delim, out_delim), args.output_file)

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser("Slice a text file. Compression will be inferred from the file names extension")
    parser.add_argument("-input_file")
    parser.add_argument("-output_file")
    parser.add_argument("-delimiter", help="String or pattern to split the input file contents (and produce the output)")
    parser.add_argument("-columns", help="Which columns to keep", nargs="+")
    args = parser.parse_args()
    run(args)