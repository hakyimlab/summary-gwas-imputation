__author__ = "alvaro barbeira"

import logging

import numpy
from scipy import stats

from ...Exceptions import ReportableException
from ...Constants import EFFECT_ALLELE, NON_EFFECT_ALLELE, ZSCORE, \
    OR, BETA, BETA_SIGN, SE, PVALUE
from ...data_management import TextFileTools


def load_gwas(path, spec, order=None, force_special_handling=False, skip_until_header=None, separator=None, handle_empty_columns=False, input_pvalue_fix=None, enforce_numeric_columns=None):
    #TODO: use snp id as key, etc
    d = TextFileTools.load_dataframe(path, spec, order=order,
        force_special_handling=force_special_handling, skip_until_header=skip_until_header, separator=separator, handle_empty_columns=handle_empty_columns)
    d = _ensure_columns(d, input_pvalue_fix, enforce_numeric_columns)
    return d

def _ensure_columns(d, input_pvalue_fix=1e-50, enforce_numeric_columns=None):
    if d.shape[0] == 0:
        if OR in d: d[BETA] = None
        if BETA_SIGN in d: d[BETA_SIGN] = None
        d[ZSCORE] = None
        return d

    if not EFFECT_ALLELE in d or not NON_EFFECT_ALLELE in d:
        logging.warning("No allele columns in GWAS! I hope you know what you are doing.")
    else:
        d[EFFECT_ALLELE] = d[EFFECT_ALLELE].str.upper()
        d[NON_EFFECT_ALLELE] = d[NON_EFFECT_ALLELE].str.upper()

    if enforce_numeric_columns:
        d = _enforce_numeric_columns(d)

    if OR in d:
        logging.log(9, "Calculating beta from odds ratio")
        beta = _or_to_beta(d[OR])
        d[BETA] = beta

    if BETA_SIGN in d:
        b = d[BETA_SIGN]
        b = b.apply(lambda x: 1.0 if x == "+" else -1.0)
        d[BETA_SIGN] = b

    _ensure_z(d, input_pvalue_fix)

    d[ZSCORE] = numpy.array(d[ZSCORE], dtype=numpy.float32)
    if BETA in d and not SE in d:
        d[SE] = d[BETA]/d[ZSCORE]

    if not PVALUE in d:
        d[PVALUE] = 2*stats.norm.sf(numpy.absolute(d[ZSCORE]))
    return d

_numeric_columns = [BETA, OR, SE, PVALUE, ZSCORE]
def _enforce_numeric_columns(d):
    logging.log(8, "Enforcing numeric columns")
    for column in _numeric_columns:
        if column in d:
            a = d[column]
            if a.dtype == numpy.object:
                a = [str(x) for x in a]
                a = [TextFileTools.sanitize_component(x) for x in a]
            d = d.assign(**{column:numpy.array(a, dtype=numpy.float64)})
    return d

def _ensure_z(d, input_pvalue_fix):
    if ZSCORE in d:
        logging.log(9, "Using declared zscore")
        return d

    z = None

    if PVALUE in d:
        logging.log(9, "Calculating zscore from pvalue")
        z = _z_from_p(d, input_pvalue_fix)
    elif SE in d and BETA in d:
        logging.info("Calculating zscore from se and beta")
        z = d[BETA] / d[SE]

    if z is None: raise ReportableException("Couldn't get zscore from GWAS")
    d[ZSCORE] = z
    return d

def _z_from_p(d, input_pvalue_fix):
    p = d[PVALUE].values
    if numpy.any(p == 0):
        logging.warning("Encountered GWAS pvalues equal to zero. This might be caused by numerical resolution. Please consider using another scheme such as -beta- and -se- columns, or checking your input gwas for zeros.")

    s = _beta_sign(d)
    abs_z = -stats.norm.ppf(p / 2)

    if numpy.any(numpy.isinf(abs_z)) and input_pvalue_fix:
        logging.warning("Applying thresholding to divergent zscores. You can disable this behavior by using '--input_pvalue_fix 0' in the command line")
        the_min = numpy.min(p[numpy.logical_and(numpy.isfinite(abs_z),p != 0)])
        if input_pvalue_fix < the_min:
            the_min = input_pvalue_fix
        fix_z = -stats.norm.ppf(the_min / 2)
        logging.warning("Using %f to fill in divergent zscores", fix_z)
        abs_z[numpy.isinf(abs_z)] = fix_z

    z = abs_z * s
    return z

def _beta_sign(d):
    b = None

    if BETA in d:
        logging.log(9, "Acquiring sign from beta")
        b = numpy.sign(d[BETA])
    elif BETA_SIGN in d:
        logging.log(9, "Acquiring sign")
        b = d[BETA_SIGN]
        b = b.apply(lambda x: 1.0 if (x =="+" or x==1.0) else -1.0)

    if b is None: raise ReportableException("No beta sign in GWAS")
    return b

def _or_to_beta(odd):
    if numpy.any(numpy.where(odd < 0)):
        raise ReportableException("Odd Ratios include negative values.")
    if numpy.any(numpy.where(odd == 0)):
        logging.warning("Odd Ratios column holds some [0] values")
    return numpy.log(odd)