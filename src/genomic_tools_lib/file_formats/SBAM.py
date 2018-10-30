__author__ = "alvaro barbeira"

import os
from .. import Utilities

def save_study(study, folder):
    """Assumes a simple study with data in pandas dataframes"""
    phenos = study.get_available_pheno_list()
    variants = study.get_variants()
    covariates = study.get_covariates()

    for pheno in phenos:
        p = os.path.join(folder, pheno) + ".txt"
        with Utilities.open_any(p, "w") as f:
            pheno_data = study.get_phenos([pheno])
            line = "pheno {} g {}\n".format(pheno, " ".join(map(str,pheno_data[pheno])))
            f.write(line)

            for v in variants.columns.values:
                line = "geno {} g {}\n".format(v, " ".join(map(str,variants[v])))
                f.write(line)

            if covariates is not None:
                for c in covariates.columns.values:
                    line = "controlled {} g {}\n".format(c, " ".join(map(str, covariates[c])))
                    f.write(line)


