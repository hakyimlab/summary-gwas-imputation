__author__ = "alvaro barbeira"

from collections import namedtuple

import pyarrow as pa

GTExAllAssociations = namedtuple("GTExAllAssociations",
["gene", "variant_id", "rsid", "maf", "pvalue", "beta", "se"])

pyarrow_schema = pa.schema([pa.field("gene", pa.string()), pa.field("variant_id", pa.string()), pa.field("rsid", pa.string()),
    pa.field("maf", pa.float32()), pa.field("pvalue", pa.float64()), pa.field("beta", pa.float32()), pa.field("se", pa.float32())])