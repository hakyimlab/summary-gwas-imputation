__author__ = "alvaro barbeira"

from genomic_tools_lib.individual_data.Genotype import _get_variants_metadata

class _Study():
    def get_genotype(self): raise RuntimeError("Not implemented")
    def get_individuals(self): raise RuntimeError("Not implemented")
    def get_phenotype(self): raise RuntimeError("Not implemented")

    def get_variants_metadata(self, variants=None): raise RuntimeError("Not implemented")
    def get_variants(self, variants=None, to_pandas=None, omit_individuals=False): raise RuntimeError("Not implemented")
    def get_available_variant_list(self): raise RuntimeError("Not implemented")
    def get_available_pheno_list(self): raise RuntimeError("Not implemented")
    def get_phenos(self, phenos=None, to_pandas=None): raise RuntimeError("Not implemented")
    def get_available_covariate_list(self): raise RuntimeError("Not implemented")
    def get_covariates(self, covariates=None, to_pandas=None): raise RuntimeError("Not implemented")

class Study(_Study):
    """
genotype is of class Genotype, phenotype is a pandas dataframe with at least one column, individuals is a list of individual ids.
    """
    def __init__(self, genotype, phenotype, individuals, covariates=None):
        self.genotype = genotype
        self.phenotype = phenotype
        self.individuals = individuals
        self.covariates = covariates

    def get_individuals(self): return self.individuals
    def get_phenotype(self): return self.phenotype
    def get_phenos(self, phenos=None, to_pandas=True): return _get(self.phenotype, phenos, to_pandas)
    def get_available_pheno_list(self): return _get_list(self.phenotype)
    def get_covariates(self, covariates=None, to_pandas=True): return  _get(self.covariates, covariates, to_pandas)
    def get_available_covariate_list(self): return _get_list(self.covariates)
    def get_genotype(self): return self.genotype
    def get_variants_metadata(self, variants=None): return _get_variants_metadata(self.genotype.metadata, variants)
    def get_variants(self, variants=None, to_pandas=True, omit_individuals=False, specific_individuals=None): return self.genotype.get_variants(variants, to_pandas, specific_individuals)

def _get(d, names=None, to_pandas=True):
    if d is None: return None
    r = d
    if names:
        r = r[names]
    if not to_pandas:
        r = {name:r[name].values for name in r.columns.values}
    return r

def _get_list(d):
    if d is None: return None
    l = [x for x  in d.columns.values if x != "individual"]
    return l
