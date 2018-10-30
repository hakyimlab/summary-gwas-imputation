__author__ = "alvaro barbeira"

import pandas
from . import KeyedDataSource

def load_gtex_variant_to_rsid(path):
    return KeyedDataSource.load_data(path, "variant_id", "rs_id_dbSNP150_GRCh38p7", should_skip=(lambda k,v: v == "."))

def load_gtex_variants(path, frequency_filter=None):
    if not frequency_filter:
        v = pandas.read_table(path, usecols=["variant_id"])
    else:
        raise RuntimeError("Frequency filter not implemented here")
    return {x for x in v}
