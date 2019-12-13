
This repository contains tools to harmonize GWAS summary statistics to a given reference.
The main application is harmonization of a public gwas' variants to those in the GTEx study,
and allow imputation of summary satistics for missing variants.

You can also find tools to run colocalization (COLOC or Enloc) on them.

# Harmonization and Imputation Overview

The first step consists of compiling a dataset with variants' definitions and genotypes in a given sample.
A script in this repository will convert them to an internal format using Apache Parquet for fast querying.

The second step is harmonizing each GWAS to the reference data set, optionally including a liftover conversion of chromosomal coordinates.

The third step is imputing missing variants' summary statistics. Because of the long computation time, this step can be split across multiple executions
to decrease total running time.

Th last step is a simple postprocessing that collects the harmonized GWAS, the imputed variants, and produces a final product.

# PRerequisites

# Other tools

Colocalization analysis via COLOC is supported on harmonized or imputed GWAS.
We include a modified copy of COLOC for ease of use (original source at [github](https://github.com/chr1swallace/coloc))
We include a modified copy of CTIMP for ease of use (original source at [github](https://github.com/yiminghu/CTIMP))
See their respective licenses at their repositories, or here in `3rd_licenses`)

# Final comments

Please see the [wiki](https://github.com/hakyimlab/summary-gwas-imputation/wiki) for extensive, in-detail documentation.
