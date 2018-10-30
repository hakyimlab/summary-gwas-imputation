__author__ = "alvaro barbeira"

import numpy
import pandas
import logging
import re

from .. import Utilities
from ..individual_data import Genotype

def _mean(x):
    return float(numpy.mean(numpy.array([_x for _x in x if _x != "NA"], dtype=numpy.float32)))

def _impute_to_mean(x):
    m = _mean(x)
    return [float(x_) if x_ != "NA" else m for x_ in x]

def dosage_generator(path, snp_key, dosage_conversion=None, dosage_filter=None, do_none=False):
    """Metadata is a tuple formatted as Genotype.MetadataTFE"""
    logging.log(9, "Opening dosage from %s", path)

    #TODO maybe add support for chromosome X
    r = re.compile("(\d+)")

    dosage_ = (lambda x: list(map(float, x))) if not dosage_conversion else dosage_conversion
    with Utilities.open_any(path) as f_:
        header_ = f_.readline().decode()
        individual_ids = header_.strip().split()[1:]
        for i,line in enumerate(f_):
            comps = line.decode().strip().split()
            variant = comps[0]

            variant_comps = variant.strip().split("_")
            chromosome = variant_comps[0]
            if "chr" in chromosome: chromosome = chromosome.split("chr")[1]

            if not r.search(chromosome): #TODO: think about this
                continue

            chromosome = int(chromosome)

            position = int(variant_comps[1])
            allele_non_effect = variant_comps[2]
            allele_effect = variant_comps[3]

            rsid = snp_key[variant] if (snp_key and variant in snp_key) else (numpy.nan if not do_none else None)

            dosage = dosage_(comps[1:])
            freq = numpy.mean(dosage) / 2
            metadata = (chromosome, position, variant, allele_non_effect, allele_effect, freq, rsid)

            if dosage_filter and dosage_filter(dosage, metadata, individual_ids):
                continue

            yield dosage, metadata, individual_ids

########################################################################################################################
def _load_genotype_file(path, snp_key, dosage_conversion, dosage_filter):
    _dosages = []
    _metadata = []
    _individual_ids = None
    for dosage, metadata, ids in dosage_generator(path, snp_key, dosage_conversion, dosage_filter):
        _dosages.append(dosage)
        _metadata.append(metadata)
        _individual_ids = ids
    return _dosages, _metadata, _individual_ids

def load_genotype_file(path, snp_key, dosage_conversion=None, filter=None):
    dosages_, metadata_, ids_ = _load_genotype_file(path, snp_key, dosage_conversion, filter)
    metadata_ = Genotype._metadata_from_raw(metadata_, extra_columns=["rsid"])
    return Genotype.Genotype(dosages_, metadata_), ids_

#######################################################################################################################

def _load_genotype_file_by_chromosome(path, snp_key, dosage_conversion, dosage_filter):
    _dosages = []
    _metadata = []
    _individual_ids = None
    last_chr=None
    for dosage, metadata, ids in dosage_generator(path, snp_key, dosage_conversion, dosage_filter):
        _chr=metadata[Genotype.MetadataTF.CHROMOSOME]

        if not last_chr:
            last_chr = _chr

        if last_chr != _chr:
            yield _dosages, _metadata, _individual_ids
            last_chr = _chr
            _dosages, _metadata, _individual_ids = [], [], []

        _dosages.append(dosage)
        _metadata.append(metadata)
        _individual_ids = ids

    yield _dosages, _metadata, _individual_ids

def load_genotype_file_by_chromosome(path, snp_key, dosage_conversion=None, filter=None):
    for _dosages, _metadata, _ids in _load_genotype_file_by_chromosome(path, snp_key, dosage_conversion, filter):
        _metadata = Genotype._metadata_from_raw(_metadata, extra_columns=["rsid"])
        yield Genotype.Genotype(_dosages, _metadata), _ids

########################################################################################################################

def load_genotype_folder(folder, pattern, snp_key):
    """Deprecated"""
    paths = sorted(Utilities.folder_contents(folder, pattern))

    dosages_, metadata_, individual_ids_ = [], [], []
    for path in paths:
        for dosage, metadata, ids in dosage_generator(path, snp_key):
            dosages_.append(dosage)
            metadata_.append(metadata)
            individual_ids_ = ids

    metadata_ = Genotype._metadata_from_raw(metadata_)
    return Genotype.Genotype(dosages_, metadata_), individual_ids_

def load_variable_file(expression_path):
    variables = pandas.read_table(expression_path)
    key = variables.columns.values[0]
    variables = variables.rename(columns={key:"name"})
    variables = variables.set_index("name").transpose().reset_index().rename(columns={"index": "individual"})
    columns = list(variables.columns.values)
    columns.remove("individual")
    columns = ["individual"] + columns
    return variables[columns]



