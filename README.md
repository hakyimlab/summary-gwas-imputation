
This repository contains tools to harmonize GWAS summary statistics to a given reference.
The main application is harmonization of a public gwas' variants to those in the GTEx study,
and allow imputation of summary statistics for missing variants.

You can also find tools to run colocalization (COLOC or Enloc) on them.

# Harmonization and Imputation Overview

The first step consists of compiling a dataset with variants' definitions and genotypes in a given sample.
A script in this repository will convert them to an internal format using Apache Parquet for fast querying.

The second step is harmonizing each GWAS to the reference data set, optionally including a liftover conversion of chromosomal coordinates.

The third step is imputing missing variants' summary statistics. Because of the long computation time, this step can be split across multiple executions
to decrease total running time.

Th last step is a simple postprocessing that collects the harmonized GWAS, the imputed variants, and produces a final product.

# Prerequisites

The basic requirements for running GWAS summary-imputation are python>=3.5 with the following packages:
 
- pandas=0.25.3
- scipy=1.4.1
- numpy=1.18.1
- bgen_reader=3.0.2
- cyvcf2=0.20.0
- pyliftover=0.4
- pyarrow=0.11.1

A quick-and-dirty solution to install these requirements is using [Miniconda](https://www.anaconda.com/open-source) and the file `src/conda_env.yaml`
in this repository to create a working environment.

```bash
conda env create -f /path/to/this/repo/src/conda_env.yaml
conda activate imlab
```

# Other tools

Colocalization analysis via COLOC is supported on harmonized or imputed GWAS.
We include a modified copy of COLOC for ease of use (original source at [github](https://github.com/chr1swallace/coloc))
We include a modified copy of CTIMP for ease of use (original source at [github](https://github.com/yiminghu/CTIMP))
See their respective licenses at their repositories, or here in `3rd_licenses`)

# Final comments

Please see the [wiki](https://github.com/hakyimlab/summary-gwas-imputation/wiki) for extensive, in-detail documentation.
