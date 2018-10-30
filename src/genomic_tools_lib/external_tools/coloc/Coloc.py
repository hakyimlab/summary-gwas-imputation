__author__ = "alvaro barbeira"

import numpy
import pandas
import logging
import traceback

import rpy2.robjects as robjects
import rpy2.robjects.numpy2ri
rpy2.robjects.numpy2ri.activate()
coloc_r = None
K_COLUMNS = None

from ... import Utilities
from ...data_management import DataFrameStreamer

#Awful hack to avoid having to learn the python import process
def initialize(script_path=None):
    global coloc_r
    global K_COLUMNS
    if coloc_r:
        raise RuntimeError("Coloc re-initialization not allowed.")
    if not script_path:
        from rpy2.robjects.packages import importr
        coloc_r = importr('coloc').coloc_abf
        K_COLUMNS = ["gene_id", "p0", "p1", "p2", "p3", "p4"]
    else:
        r = robjects.r
        r['source'](script_path)
        coloc_r = r['coloc.abf']
        K_COLUMNS = ["gene_id", "p0", "p1", "p2", "p3", "p4", "cp0", "cp1", "cp2", "cp3", "cp4"]

########################################################################################################################

def _sanitize(d):
    if "pvalue" in d:
        d = d.loc[d.pvalue>0]
    if "frequency" in d:
        d = d.loc[d.frequency>0]
    if "maf" in d:
        d = d.loc[d.maf>0]
    return d

########################################################################################################################
def _read(gwas, cols, gwas_sample_size):
    if gwas_sample_size == "FROM_GWAS":
        cols += ["sample_size"]
    d = pandas.read_table(gwas, usecols=cols)
    logging.log(9, "sanitizing gwas")
    #d = _sanitize(d)
    gwas_sample_size = None if gwas_sample_size == "FROM_GWAS" else int(gwas_sample_size)
    return d, gwas_sample_size

def _read_gwas_pvalue(gwas, gwas_sample_size):
    d, gwas_sample_size = _read(gwas, ["panel_variant_id", "pvalue", "frequency"], gwas_sample_size)
    return {x.panel_variant_id: (x.pvalue, x.frequency, gwas_sample_size if gwas_sample_size else x.sample_size) for x in d.itertuples() if (x.pvalue>0 and x.frequency>0)}

def _read_gwas_bse(gwas, gwas_sample_size):
    d, gwas_sample_size = _read(gwas, ["panel_variant_id", "effect_size", "standard_error", "frequency"], gwas_sample_size)
    return {x.panel_variant_id: (x.effect_size, x.standard_error, x.frequency, gwas_sample_size if gwas_sample_size is not None else x.sample_size)for x in d.itertuples() if x.frequency>0}

def _read_gwas_zscore_1(gwas, gwas_sample_size):
    d, gwas_sample_size = _read(gwas, ["panel_variant_id", "zscore", "frequency"], gwas_sample_size)
    return { x.panel_variant_id: (x.zscore, 1.0, x.frequency, gwas_sample_size if gwas_sample_size is not None else x.sample_size) for x in d.itertuples() if x.frequency>0}

def read_gwas(gwas, gwas_sample_size, gwas_mode="pvalue"):
    methods = {"pvalue": _read_gwas_pvalue, "bse": _read_gwas_bse, "zscore_1":_read_gwas_zscore_1}
    method = methods.get(gwas_mode)
    if not method: raise RuntimeError("unsupported gwas mode")

    return method(gwas, gwas_sample_size)

########################################################################################################################
def eqtl_streamer(eqtl_path, keys):
    columns = ["maf", "pval_nominal", "slope", "slope_se"]
    #_skip = lambda x: x not in keys
    #return DataFrameStreamer.data_frame_streamer(eqtl_path, sanitize=True, to_numeric=columns, sentinel_column="gene_id", additional_skip_row_check=_skip)
    #
    return DataFrameStreamer.data_frame_streamer(eqtl_path, sanitize=True, to_numeric=columns, sentinel_column="gene_id")

def get_eqtl(d, eqtl_sample_size, eqtl_mode="bse"):
    logging.log(9, "sanitizing eqtl")

    ####################################################################################################################
    # This nice code can't be because of an obscure coloc bug, concerning MAF>0 check.
    # Coloc thinks itself smart and able to handle MAF==0 but it fails miserably.
    # if eqtl_mode == "pvalue": l = lambda x: (x.pval_nominal, x.maf, eqtl_sample_size)
    # elif eqtl_mode == "bse": l = lambda x: (x.slope, x.slope_se, x.maf, eqtl_sample_size)
    # elif eqtl_mode == "zscore_1": l = lambda x: (x.slope/x.slope_se, 1.0, x.maf, eqtl_sample_size)
    # else: raise RuntimeError("unsupported eqtl mode")
    #return {x.variant_id: l(x) for x in d.itertuples()}

    ####################################################################################################################
    # you suck bigtime, coloc. May a ciphrang feast on you on the outside.
    if eqtl_mode == "pvalue":
        r = {x.variant_id:(x.pval_nominal, x.maf, eqtl_sample_size) for x in d.itertuples() if (x.pval_nominal>0 and x.maf>0)}
    elif eqtl_mode == "bse":
        r = {x.variant_id:(x.slope, x.slope_se, x.maf, eqtl_sample_size) for x in d.itertuples() if x.maf>0}
    elif eqtl_mode == "zscore_1":
        r = {x.variant_id: (x.slope/x.slope_se, 1.0, x.maf, eqtl_sample_size) for x in d.itertuples() if x.maf>0}
    else:
        raise RuntimeError("unsupported eqtl mode")
    return r

########################################################################################################################

def coloc_on_gwas_eqtl(gene, gwas, eqtl, gwas_mode, eqtl_mode, p1=1e-4, p2=1e-4, p12=1e-5):
    g = {k: gwas[k] for k in eqtl if k in gwas}
    keys = sorted(g.keys())
    _gwas = _to_coloc_data(g, gwas_mode, keys)
    _eqtl = _to_coloc_data(eqtl, eqtl_mode, keys)
    return _coloc(gene, _gwas, _eqtl, p1, p2, p12)

def _convert(d, mode, keys):
    if mode == "pvalue":
        #converted = pandas.DataFrame([d[k] for k in keys], columns=["pvalue", "frequency", "sample_size"])
        converted = pandas.DataFrame([(k,)+d[k] for k in keys], columns=["id", "pvalue", "frequency", "sample_size"])
    elif mode == "bse" or mode == "zscore_1":
        # zscore_1 has a different meaning but all processing difference is upstream
        #converted = pandas.DataFrame([d[k] for k in keys], columns=["beta", "se", "frequency", "sample_size"])
        converted = pandas.DataFrame([(k,)+d[k] for k in keys], columns=["id", "beta", "se", "frequency", "sample_size"])
    else:
        raise RuntimeError("unsupported mode")
    return converted

def _to_coloc_data(d, mode, keys, type="quant"):
    converted = _convert(d, mode, keys)
    return  d_(converted, None, type)

########################################################################################################################
def _s(d, N):
    if "sample_size" in d:
        return d["sample_size"].values
    elif N:
        return N
    else:
        raise RuntimeError("Need sample size for p-value-based coloc")

def d_p_(d, N=None, type="quant"):
    p = robjects.FloatVector(d["pvalue"].values)
    f = robjects.FloatVector(d["frequency"].values)
    s = _s(d, N)
    return robjects.ListVector({"pvalues": p, "MAF": f, "N": s, "type": type})

def d_bse_(d, N, type="quant"):
    v = robjects.FloatVector(numpy.array(d["se"]) ** 2)
    f = robjects.FloatVector(d["frequency"].values)
    b = robjects.FloatVector(d["beta"].values)
    s = _s(d, N)
    return robjects.ListVector({"beta": b, "varbeta": v, "MAF": f, "N": s, "type": type})

def d_(d, N, type="quant"):
    if "pvalue" in d:
        return d_p_(d, N, type)
    elif "beta" in d and "se" in d:
        return d_bse_(d, N, type)
    else:
        raise RuntimeError("Data not supported for coloc")

########################################################################################################################

def coloc_from_dataframe(gene, df1, df2, N1=None, N2=None, p1=1e-4, p2=1e-4, p12=1e-5):
    d1 = d_(df1, N1)
    d2 = d_(df2, N2)
    return _coloc(gene, d1, d2, p1, p2, p12)

def _coloc(gene, d1, d2, p1, p2, p12):
    try:
        c = coloc_r(dataset1=d1, dataset2=d2, p1=p1, p2=p2, p12=p12)
        r = c.rx('summary')[0]
        r = tuple([gene]+list(r[1:]))
    except Exception as ex:
        logging.info("Exception running coloc:\n%s", traceback.format_exc())
        r = tuple([None]*len(K_COLUMNS))
    return r

def results_to_dataframe(r):
    return Utilities.to_dataframe(r, K_COLUMNS)
