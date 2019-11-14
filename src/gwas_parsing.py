#!/usr/bin/env python
__author__ = "alvaro barbeira"
import logging
import os
import math
from timeit import default_timer as timer

import numpy
import pandas
import pyliftover

from genomic_tools_lib import Logging, Utilities
from genomic_tools_lib.miscellaneous import PandasHelpers, Genomics
from genomic_tools_lib.data_management import TextFileTools
from genomic_tools_lib.file_formats.gwas import GWAS, Utilities as GWASUtilities

def try_convert(x, t):
    c=True
    try:
        x = t(x)
    except:
        c = False
    return c

def inferred(x):
    if try_convert(x, int):
        return int(x)
    elif try_convert(x, float):
        return float(x, float)

    return x

def insert_value(d, value_spec):
    x = inferred(value_spec[1])
    d[value_spec[0]] = x
    return d

def pre_process_gwas(args, d):
    if args.split_column:
        for s in args.split_column:
            d = PandasHelpers.split_column(d, s)
        if "position" in d:
            d = d.assign(position = d.position.astype(int))

    if args.insert_value:
        for spec in args.insert_value:
            d = insert_value(d, spec)

    # Some GWAs have NA's in fre
    if "frequency" in d:
        d["frequency"] = Genomics.to_number(d.frequency)

    if "n_controls" in d:
        if "n_cases" in d:
            logging.info("Adding up to sample size")
            d["sample_size"] = d.n_cases + d.n_controls
        elif "sample_size" in d:
            logging.info("difference to cases")
            d["n_cases"] = d.sample_size - d.n_controls

    return d

def metadata_white_list(black_list_path, column, variants):
    w = {x for x in variants}
    if black_list_path:
        b = TextFileTools.load_column(black_list_path, column, unique_entries=True, white_list=w)
        w = {x for x in w if not x in b}
    return w

def fill_coords(args, d):
    logging.info("Loading SNP metadata whitelist")
    w = metadata_white_list(args.snp_info_blacklist, "name", d.variant_id)

    logging.info("Loading SNP specification")
    s = TextFileTools.load_dataframe(args.fill_from_snp_info, keys=w, key_column_name="name").rename(columns={"start":"position"})

    d = d.merge(s, left_on="variant_id", right_on="name", how="left")
    logging.info("%d variants after filling coordinates", d.shape[0])
    return d

def _lift(l, chromosome, position):
    # NA is important, instead of None or NaN, so that integer positions are not converted to floats by pandas. Yuck!
    _new_chromosome = "NA"
    _new_position = "NA"
    try:
        p = int(position)
        l_ = l.convert_coordinate(chromosome, p)
        if l_:
            if len(l_) > 1:
                logging.warning("Liftover with more than one candidate: %s", t.variant_id)
            _new_chromosome = l_[0][0]
            _new_position = int(l_[0][1])
    except:
        pass
    return _new_chromosome, _new_position

def liftover(args, d):
    logging.info("Performing liftover")
    l = pyliftover.LiftOver(args.liftover)
    new_position = []
    new_chromosome = []
    for t in d.itertuples():
        _new_chromosome, _new_position = _lift(l, t.chromosome, t.position)

        new_chromosome.append(_new_chromosome)
        new_position.append(_new_position)

    d = d.assign(chromosome=new_chromosome)
    d = d.assign(position=new_position)

    logging.info("%d variants after liftover", d.shape[0])
    return d

def _get_panel_metadata(path, index):
    m=[]
    for i, line in Utilities.iterate_file(path):
        if i==0: continue
        comps = line.strip().split()
        variant = comps[0]
        chr = comps[1]
        pos = int(comps[2])
        non_effect = comps[3]
        effect = comps[4]
        frequency = float(comps[5])
        if chr in index and pos in index[chr]:
            m.append((variant, chr, pos, non_effect, effect, frequency))
    return m

def _get_metadata(path, index):
    m=[]
    for i, line in Utilities.iterate_file(path):
        if i==0: continue
        comps = line.strip().split()
        chr = "chr"+comps[0]
        pos = int(comps[1])
        variant = comps[2]
        non_effect = comps[3]
        effect = comps[4]
        frequency = float(comps[5]) if comps[5] != "NA" else numpy.nan
        if chr in index and pos in index[chr]:
            m.append((variant, chr, pos, non_effect, effect, frequency))
    return m

def get_panel_variants(args, d, keep_rsids=False):
    logging.info("Creating index to attach reference ids")
    index = {}
    for t in d.itertuples():
        if not t.chromosome in index: index[t.chromosome] = set()
        index[t.chromosome].add(t.position)

    logging.info("Acquiring reference metadata")
    _snp = args.snp_reference_metadata
    m = None
    if len(_snp) == 1:
        m = _get_panel_metadata(_snp[0], index)
    elif len(_snp) == 2:
        if _snp[1] == "METADATA":
            m = _get_metadata(_snp[0], index)
    m = pandas.DataFrame(m, columns=["panel_variant_id", "chromosome", "position", "non_effect_allele", "effect_allele", "frequency"])
    return m

def filled_frequency(d, m):
    has_freq = "frequency" in d

    t = {x.panel_variant_id:x.frequency for x in m.itertuples()}
    f = []
    for e in d.itertuples():
        _f = e.frequency if has_freq else None
        if _f is None or numpy.isnan(_f):
            if e.panel_variant_id in t:
                _f = t[e.panel_variant_id]
        f.append(_f)
    return f

def _ensure_uniqueness(d):
    top_ = d[["panel_variant_id", "zscore", "variant_id"]].assign(absz = numpy.abs(d.zscore)).drop(["zscore"], axis=1)
    top_ = top_.groupby(["panel_variant_id"]).\
        apply(lambda x: x.sort_values(by="absz", ascending=False)).\
        reset_index(drop=True).groupby(["panel_variant_id"]).\
        head(1)
    top_variants = {x for x in top_.variant_id}
    exclude = {x for x in d.variant_id if not x in top_variants}
    d = d[~d.variant_id.isin(exclude)]
    return d

def ensure_uniqueness(d):
    #ugly code to avoid pandas inefficiency in grouping operations
    # key = []
    # for e in d.itertuples():
    #     _key = "NA" if not e.panel_variant_id else e.panel_variant_id
    #     key.append(_key)
    # d["kk"] = key

    v = {}
    _v = set()
    for e in d.itertuples():
        _i = e.panel_variant_id
        if _i in v:
            _v.add(_i)
        else:
            v[_i] = []
        v[_i].append((_i, e.Index, math.fabs(e.zscore)))

    remove=set()
    for id in _v:
        entries = v[id]
        largest = entries[0]

        for i in range(1, len(entries)):
            entry = entries[i]
            if entry[2] > largest[2]:
                largest = entry

        for e in entries:
            if e[1] != largest[1]:
                remove.add(e[1])

    d = d.drop(remove)
    return d

def fill_from_metadata(args, d):
    m = get_panel_variants(args, d)
    if "panel_variant_id" in d: d= d.drop(["panel_variant_id"])

    logging.info("alligning alleles")

    d = Genomics.match(d, m)

    if not args.keep_all_original_entries:
        d = d.loc[~d.panel_variant_id.isna()]
        logging.info("%d variants after restricting to reference variants", d.shape[0])

        logging.info("Ensuring variant uniqueness")
        d = ensure_uniqueness(d)
        logging.info("%d variants after ensuring uniqueness", d.shape[0])

    logging.info("Checking for missing frequency entries")
    d["frequency"] = filled_frequency(d, m)

    return d

def clean_up(d):
    d = d.assign(sample_size=[int(x) if not math.isnan(x) else "NA" for x in d.sample_size])
    if "chromosome" in d.columns.values and "position" in d.columns.values:
        d = Genomics.sort(d)
    return d

def run(args):
    if os.path.exists(args.output):
        logging.info("output path %s exists. Nope.", args.output)
        return

    start = timer()
    logging.info("Parsing input GWAS")
    d = GWAS.load_gwas(args.gwas_file, args.output_column_map,
            force_special_handling=args.force_special_handling, skip_until_header=args.skip_until_header,
            separator=args.separator, handle_empty_columns=args.handle_empty_columns, input_pvalue_fix=args.input_pvalue_fix,
            enforce_numeric_columns=args.enforce_numeric_columns)
    logging.info("loaded %d variants", d.shape[0])

    d = pre_process_gwas(args, d)

    if args.fill_from_snp_info:
        d = fill_coords(args, d)

    if args.chromosome_format:
        d = d.assign(chromosome = Genomics.to_int(d.chromosome))
        d = d.assign(chromosome = ["chr{}".format(x) for x in d.chromosome])

    if args.liftover:
        d = liftover(args, d)

    if args.snp_reference_metadata:
        d = fill_from_metadata(args, d)

    if args.output_order:
        order = args.output_order
        for c in order:
            if not c in d:
                d= d.assign(**{c:numpy.nan})
        d = d[order]

    d = clean_up(d)

    logging.info("Saving...")
    Utilities.save_dataframe(d, args.output, fill_na=True)
    end = timer()
    logging.info("Finished converting GWAS in %s seconds", str(end-start))

if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser("parse a GWAS into a standard format")
    parser.add_argument("-gwas_file", help="GWAS summary statistics file")
    parser.add_argument("-liftover", help = "File with liftover chain")
    parser.add_argument("-snp_reference_metadata", help="File with reference (GTEx) snp metadata", nargs="+")
    parser.add_argument("-snp_info_blacklist", help="ignore snps in this list")
    parser.add_argument("-split_column", help="Specify multiple key-value pairs to specify format conversion", nargs="+", action="append")
    parser.add_argument("--chromosome_format", help="Convert chromosome column to -chr{}- ", action="store_true")
    parser.add_argument("-output_column_map", help="Specify multiple key-value pairs to specify format conversion", nargs=2, action="append")
    parser.add_argument("--insert_value", help="Create a column with a specific value", nargs=2, action="append")
    parser.add_argument("-output_order", help="Specify output order", nargs='+')
    parser.add_argument("-output", help="Where the output should go")
    parser.add_argument("--keep_all_original_entries", action="store_true")
    parser.add_argument("-verbosity", help="Log verbosity level. 1 is everything being logged. 10 is only high level messages, above 10 will hardly log anything", default = "10")
    GWASUtilities.add_gwas_arguments_to_parser(parser)
    args = parser.parse_args()

    Logging.configure_logging(int(args.verbosity))

    run(args)
