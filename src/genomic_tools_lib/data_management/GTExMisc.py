__author__ = "alvaro barbeira"

import pandas
from . import KeyedDataSource

def load_gtex_variant_to_rsid(path, rsid_column = "rs_id_dbSNP150_GRCh38p7"):
    return KeyedDataSource.load_data(path, "variant_id", rsid_column, should_skip=(lambda k,v: v == "."))

def load_gtex_variants(path, frequency_filter=None):
    if not frequency_filter:
        v = pandas.read_table(path, usecols=["variant_id"])
    else:
        raise RuntimeError("Frequency filter not implemented here")
    return {x for x in v}
