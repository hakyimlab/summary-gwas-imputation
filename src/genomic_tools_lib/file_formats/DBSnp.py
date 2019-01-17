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

def column_index(recode_observed=False):
    index = {x:i for i,x in enumerate(DBSNP._fields)}
    if recode_observed:
        _o = index["observed"]
        k = {k:(v+1 if v > index["observed"] else v) for k,v in index.items()}
        del index["observed"]
        index["reference"] = _o
        index["alternate"] = _o+1
    return index

_variant_index = column_index()["name"]
_chr_index_ = column_index()["chromosome"]

def _report(variant, reported, repeat_message=""):
    if variant in reported: return
    logging.log(7, "Skipping %s %s", repeat_message, variant)
    reported.add(variant)

def _shifted(comps, index, field, shift=1):
    _i = index[field]
    f = int(comps[_i])
    return str(f + shift)

def _generate(input_file, keep_zero_based=False, recode_observed=True):
    index = column_index()

    for i, line in Utilities.iterate_file(input_file):
        comps = line.strip().split()
        if not keep_zero_based:

            comps[index["start"]] = _shifted(comps, index, "start")
            if comps[index["variant_class"]] == "deletion":
                comps[index["end"]] = _shifted(comps, index, "end")

            if comps[index["variant_class"]] == "insertion":
                comps[index["end"]] = _shifted(comps, index, "end", 2)

        if recode_observed:
            _i = index["observed"]
            _o = comps[_i]

            if "/" in _o:
                _o = _o.split("/")
                if len(_o) > 2:
                    logging.log(9, "Skipping unknown observed: %s", comps[_i])
                    #continue
                if comps[index["variant_class"]] == "deletion":
                    _r = _o[1]
                    _a = _o[0]
                else:
                    _r = _o[0]
                    _a = _o[1]
                del comps[_i]
                comps.insert(_i, _r)
                comps.insert(_i+1, _a)
            else:
                _r="error"
                _a="error"
        yield i, comps

def _indexes(fields, recode_observed=False):
    index = column_index(recode_observed)
    indexes = [index[k] for k in fields]
    return indexes

def generate(input_file, fields, keep_zero_based=False, black_list=None, only_valid_chromosomes=True, recode_observed=False):
    indexes = _indexes(fields, recode_observed)
    for i, comps in _generate(input_file, keep_zero_based, recode_observed):
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

def generate_skips(input_file, fields, keep_zero_based=False, recode_observed=False):
    indexes = _indexes(fields, recode_observed)

    found = set()
    repeated = set()

    for i, comps in _generate(input_file, keep_zero_based, recode_observed):
        if not is_valid_chr(comps):
            continue

        variant = comps[_variant_index]
        if variant in found:
            repeated.add(variant)
            continue

        found.add(variant)

    for i, comps in _generate(input_file, keep_zero_based, recode_observed):
        variant = comps[_variant_index]
        if variant in repeated:
            yield i, [comps[_i] for _i in indexes] + ["repeated"]
            continue

        if not is_valid_chr(comps):
            yield i, [comps[_i] for _i in indexes] + ["discarded_chromosome"]




