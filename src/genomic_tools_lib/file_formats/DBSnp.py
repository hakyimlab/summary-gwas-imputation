__author__ = "alvaro barbeira"

import re
import logging
from collections import namedtuple

from .. import Utilities

DBSNP = namedtuple("SBSNP",
["binary_index", "chromosome", "start","end", "name",  "score", "strand", "ref_ncbi", "ref_ucsc", "observed",
 "mol_type", "variant_class", "valid", "average_heterozigosity", "average_heterozigosity_se", "functional_category",
 "location_type", "quality_weight", "exceptions", "submitter_count", "submitters", "alleles_with_freq",
 "alleles", "allele_chromosome_count", "allele_frequencies", "bitfields"])

def column_index():
    return {x:i for i,x in enumerate(DBSNP._fields)}

_variant_index = column_index()["name"]
_chr_index_ = column_index()["chromosome"]

def _report(variant, reported, repeat_message=""):
    if variant in reported: return
    logging.log(7, "Skipping %s %s", repeat_message, variant)
    reported.add(variant)

def _generate(input_file, keep_zero_based=False):
    index = column_index()

    for i, line in Utilities.iterate_file(input_file):
        comps = line.strip().split()
        if not keep_zero_based:
            _i = index["start"]
            start = int(comps[_i])
            start = str(start+1)
            comps[_i] = start

        yield i, comps

def generate(input_file, fields, keep_zero_based=False, black_list=None, only_valid_chromosomes=True):
    index = column_index()
    indexes = [index[k] for k in fields]
    for i, comps in _generate(input_file, keep_zero_based):
        if black_list and comps[_variant_index] in black_list:
            continue
        if only_valid_chromosomes and not is_valid_chr(comps):
            continue
        yield i, [comps[_i] for _i in indexes]

r_ = re.compile(r"chr[\d]{1,2}$")
rX_ = re.compile(r"chr[a-zA-Z]{1}$")

def is_valid_chr(comps):
    chromosome = comps[_chr_index_]
    return r_.search(chromosome) or rX_.search(chromosome)

def generate_skips(input_file, fields, keep_zero_based=False):
    index = column_index()
    indexes = [index[k] for k in fields]

    found = set()
    repeated = set()

    for i, comps in _generate(input_file, keep_zero_based):
        if not is_valid_chr(comps):
            continue

        variant = comps[_variant_index]
        if variant in found:
            repeated.add(variant)
            continue

        found.add(variant)

    for i, comps in _generate(input_file, keep_zero_based):
        variant = comps[_variant_index]
        if variant in repeated:
            yield i, [comps[_i] for _i in indexes] + ["repeated"]
            continue

        if not is_valid_chr(comps):
            yield i, [comps[_i] for _i in indexes] + ["discarded_chromosome"]




