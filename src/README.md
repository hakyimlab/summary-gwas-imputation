# Table of Contents

1. [Data Formats](#data-formats)
    1. [Variant IDs](#variant-ids)
    1. [Model Training Format](#model-training-format)
    1. [Parquet Genotype](#parquet-genotype)
    1. [Parquet Phenotype](#parquet-phenotype)
1. [Scripts](#scripts)
    1. [covariance_for_model.py](#covariance_for_modelpy)
    1. [generate_qq_plots.py](#generate_qq_plotspy)
    1. [gwas_parsing.py](#gwas_parsingpy)
    1. [model_training.py](#model_trainingpy)
    1. [model_training_genotype_to_parquet.py](#model_training_genotype_to_parquetpy)
    1. [post_process_model_training.py](#post_process_model_trainingpy)
    1. [regress_out_covariates.py](#regress_out_covariatespy)
    1. [run_susier_2.py](#run_susier_2py)
    1. [run_tensorqtl.py](#run_tensorqtlpy)
    1. [validate_model.py](#validate_modelpy)


# Data Formats

Here are a few data formats that are used frequently in this library:

## Variant IDs

Many of the scripts here are agnostic to the format of variant ID, but we have used one particular convention for the `GTEx` and `ImageXcan` projects. 

For `GTEx`, that convention is `chr<chr>_<position>_<non-effect-allele>_<effect-allele>_b38`. 

For `ImageXcan`, that convention is `chr<chr>_<position>_<non-effect-allele>_<effect-allele>`.


## Model Training Format

"Model training format" is a g-zipped text file with tabs separating columns and a variant on each row, and an individual on each column. The first column must be called `varID` and list the variant IDs. The other column headers are the individual IDs. The values should be dosages of the effect allele. 

```
varID   ID1 ID2 ID3 ... 
chr1_13526_C_T_b38  0   1   0   ...
```

Accompanying the variant file should be a SNP metadata (or annotations) file. I'm not too sure about the naming conventions for the columns. Here is what the top of my SNP annotation file looks like: 

```
chromosome	pos	varID	ref_vcf	alt_vcf	R2	MAF	rsid	rsid_dbSNP150
chr22	16050822	chr22_16050822_A_G	A	G	NA	NA	NA	rs12172168
chr22	16051249	chr22_16051249_C_T	C	T	NA	NA	NA	rs62224609
chr22	16052080	chr22_16052080_A_G	A	G	NA	NA	NA	rs4965031
chr22	16052167	chr22_16052167_AAAAC_A	AAAAC	A	NA	NA	NA	rs375684679
chr22	16052463	chr22_16052463_C_T	C	T	NA	NA	NA	rs587646183
chr22	16052684	chr22_16052684_C_A	C	A	NA	NA	NA	rs139918843
chr22	16052962	chr22_16052962_T_C	T	C	NA	NA	NA	rs376238049
chr22	16052986	chr22_16052986_A_C	A	C	NA	NA	NA	rs200777521
chr22	16053444	chr22_16053444_T_A	T	A	NA	NA	NA	rs80167676
```
## Parquet Genotype

The [Apache Parquet](https://parquet.apache.org/) format allows us to load variant data quickly (especially small groups of SNPs). There are two complementary files, one with the individual-level variant data, and one with metadata for the variants.

### Genotype File

The individual-level variant data has a simple schema: there is a column of individual IDs, named `individual`, and the other columns are named by variant IDs. See the [Variant IDs](#variant-ids) section. The `individual` column should have text individual IDs, and the variant columns should have the corresponding dosages (floats) in the order set by the `individual` column.

### Metadata File

Parquet file with the following columns: `chromosome`, `position`, `id`, `allele_0`, `allele_1`, `allele_1_frequency`, `rsid`. Here is a snapshot of the metadata file for chromosome 22 for the ImageXcan project. I don't think `allele_1_frequency` is absolutely necessary.

```IPython
In [2]: from pyarrow import parquet as pq

In [3]: fp = '/vol/bmd/meliao/data/genotype/parquet/chr-22.variants_metadata.parquet'

In [4]: df = pq.read_table(fp).to_pandas()

In [5]: df.head()
Out[5]:
   chromosome  position                      id allele_0 allele_1  allele_1_frequency         rsid
0          22  16050822      chr22_16050822_G_A        G        A            0.149176   rs12172168
1          22  16051249      chr22_16051249_T_C        T        C            0.098814   rs62224609
2          22  16052080      chr22_16052080_G_A        G        A            0.127218    rs4965031
3          22  16052167  chr22_16052167_A_AAAAC        A    AAAAC            0.438541         None
4          22  16052962      chr22_16052962_C_T        C        T            0.088033  rs376238049
```

## Parquet Phenotype

The parquet phenotype specification is a lot like the parquet genotype specification. There is one text column for individual IDs named `individual`. The other columns are the phenotype names. In the ImageXcan project, the phenotype names are `IDP-<UK Biobank Field ID #>`. 

# Scripts

Here are usage notes and command-line arguments for the scripts I've used in the ImageXcan pipeline. Some of this is legacy code, so if I don't mention a certain script, it's because I haven't worked with it. 

Mostly, the convention is that a single dash before the argument like `-this` means the argument is mandatory, and a double dash like `--this` means the argument is optional. 

## covariance_for_model.py

This is a script for creating a covariance file for a pre-existing prediction model. It takes as input genotype data in parquet format and a PredictDB prediction model. 

| Argument | Description | Example Usage |
| --- | --- | --- |
| -parquet_genotype_folder |  Folder where genotype data resides. | |
| -parquet_genotype_pattern | Filename pattern for the genotype file(s). This needs to be a regular expression. | chr-(.)_variants.parquet |
| -model_db | Path to PredictDB model. The list of weights (and corresponding SNPs) is read from this file. | |
| -output | Filepath for the output covariance file | | 
| --output_rsids | Specify this if rsids should be used rather than the normal variant ID specification | | 
| --metadata_lift | Use this to specify the path of a UCSC liftover file if necessary. Must be used with `--parquet_metadata` | | 
| --parquet_metadata | A SNP metadata file. Necessary if a liftover file is being used. | | 
| --individuals | Filepath for a list of individual IDs. Use this if only part of the genotype data should be used, i.e. some cohort of individuals. | | 
| -parsimony | Logging amount. 1 means that everything is being logged; 10 (the default) is only high-level messages | |

## generate_qq_plots.py

This script is meant to take the output of the [run_tensorqtl.py](#run_tensorqtlpy) script and plot QQ plots for all of the phenotypes. It is specific to the output directory structure of the run_tensorqtl.py script. 


| Argument | Description | Example Usage |
| --- | --- | --- |
| -out_dir | Directory to save QQ plot .png files | | 
| -in_dir | Base directory of TensorQTL GWAS results | | 
| --n_batches | Instead of running all of the phenotypes, split into this many batches. | |
| --batch | Run this batch (0-indexed) | |
| --pheno_map | If for some reason the phenotype names were assigned numbers, you can recover the pheno names with a tsv file with columns pheno_name and pheno_num | | 
| -parsimony | Logging amount. 1 means that everything is being logged; 10 (the default) is only high-level messages | |

## gwas_parsing.py

This script is explained in great detail in the wiki. The only change I've added is the `--canonical_variant_id` option, which, when specified, builds a column (with name specified after `--canonical_variant_id`) of the `chr<chromosome>_<position>_<non-effect-allele>_<effect-allele>` variant IDs.

## model_training.py

This script takes genotype data and phenotype data, as well as optional candidate SNPs, and trains ElasticNet prediction models using R's [glmnet](https://cran.r-project.org/web/packages/glmnet/index.html) package. It outputs text files of the weights table, the extra table, and the covariance matrix. Two scripts in the pipeline, [validate_model.py](#validate_modelpy) and [post_process_model_training.py](#post_process_model_trainingpy) both expect the outputs as written by this script. 

This script has a ton of arguments; I haven't used all of them. The changes I've made have likely broken some of the arguments or rendered some of them redundant, so I'll only list the ones I've used in the ImageXcan project. 

| Argument | Description | Example Usage |
| --- | --- | --- |
| -features | Either a single parquet genotype file, or a parquet genotype file that is configurable with a `chr` key. | ~/genotype/chr-{chr}_variants.parquet |
| -feautes_annotation | Either a single parquet variant metadata file, or a parquet metadata file that is configurable with a `chr` key. | ~/genotype/chr-{chr}_variants_metadata.parquet |
| -data | Phenotype data in parquet format | | 
| -output_prefix | Filename prefix. '_weights.txt.gz', '_summary.txt.gz' and '.txt.gz' will be added to the prefix to create the three output files. | |
| --mode | Either "elastic_net" (default) where glmnet is used, or "ols" where a python OLS engine is used. | | 
| --preparsed_weights | tsv with columns `gene`, `variant_id` and `weight` which creates a list of candidate SNPs for each gene. The weights are applied to the ElasticNet penalty term. Because this is expected to be finemapping output the penalty is actually `1 - weight` for all of the weights. | | 
| --preparsed_weights | tsv with columns `gene` and `variant_id` which creates a list of candidate SNPs for each gene. Specifying this option does not impose any pre-weighting of the SNPs.  | | 
| --n_train_test_folds | Number of train/test folds for model performance evaluation. Default = 5. | |
| --n_k_folds | Number of inner training folds for CV model training. Default = 10 | |
| --alpha | Specifies the alpha coefficient in the ElasticNet penalty term (which controls the $L_1 - L_2$ mixing parameter. Default is 0.5. Alpha of 0.0 is Ridge regression. | | 
| --sub_batches | Instead of running all of the phenotypes, split into this many batches. | |
| --sub_batch | Run this batch (0-indexed) | |
| -in_dir | Base directory of TensorQTL GWAS results | | 
| --n_batches | Instead of running all of the phenotypes, split into this many batches. | |
| --batch | Run this batch (0-indexed) | |
| -parsimony | Logging amount. 1 means that everything is being logged; 10 (the default) is only high-level messages | |

## model_training_genotype_to_parquet.py 

For converting a genotype file in [model training format](#model-training-format) to [parquet format](#parquet-genotype). Needs a SNP annotation file and a genotype file in model training format. 


| Argument | Description | Example Usage |
| --- | --- | --- |
| -snp_annotation_file | Filepath with SNP annotations. | | 
| -input_genotype_file | Filepath with individual-level genotype data | | 
| --only_in_key | Only keep the variants that appear in the SNP annotation file | |
| --biallelic_only | Only keep biallelic variants | |
| --filter_maf | Filters variants by minor allele frequency | | 
| --impute_to_mean | If specified, imputes missing dosages to the mean | |
| -output_prefix | Filename prefix. Will append `.variants.parquet` and `.variants_metadata.parquet` | | 
| --split_by_chromosome | Break up a big file into one file per chromosome | | 
| -parsimony | Logging amount. 1 means that everything is being logged; 10 (the default) is only high-level messages | |

## post_process_model_training.py

After running the [model_training.py](#model_trainingpy) script (and possibly after running the [validate_model.py](#validate_modelpy) script), this can be run to merge all model training output files into a PredictDB prediction model database, a covariance file, and a validation report file. 

| Argument | Description | Example Usage |
| --- | --- | --- |
| -input_prefix | The filepath prefix (also a regular expression) | ~/model_training/out/dapg_ridge_batch-(.*) | 
| -output_prefix | Filename prefix. Will have `.db`, `.txt.gz`, and `_validation.txt.gz` appended to this value | |
| -output_prefix_text | Filename prefix. The sqlite database tables will be written to text files if this is specified. | |
| -coalesce_validation | Specify if the validation files should be searched and collected | |
| --sample_info | A quick way to specify n_samples, population, tissue | 5503 YRI whole_blood |
| --skip_check | If this is not specified, some gene models will be discarded due to prediction perfomance criteria | |
| -parsimony | Logging amount. 1 means that everything is being logged; 10 (the default) is only high-level messages | |

## regress_out_covariates.py

This script takes two files formatted with the parquet phenotype format. It regresses the covariates out of the data and then writes the residuals into a new file. 

| Argument | Description | Example Usage |
| --- | --- | --- |
| -data | Phenotype data in parquet format | |
| -covariate | Covariate data in parquet format | | 
|  --inverse_normalize_data | Specify this if the phenotypes should be inverse normalized before regression | |
| -output | Filepath where the output should be written  | |
| -parsimony | Logging amount. 1 means that everything is being logged; 10 (the default) is only high-level messages | |


## run_susier_2.py

This script loads genotype and phenotype data, and it uses the [susieR package](https://stephenslab.github.io/susieR/) to compute posterior inclusion probabilities.

| Argument | Description | Example Usage |
| --- | --- | --- |
| -data | Phenotype data in parquet format | |
| -covariate | Covariate data in parquet format | | 
| -out_dir | Directory to print results to | |
| -region | Filepath for an LD Region file which defines LD-independent regions of the genome (Pickrell regions) | |
| -metadata_pattern | Parquet metadata files, formattable with a `{chr}` argument. | |
| -geno_pattern | Parquet variant files, formattable with a `{chr}` argument. | |
| --whitelist | Filepath to a text file with 'region_id' and 'gene' columns. If this is specified, only the listed region x gene pairs will be evaluated. Otherwise, all region x gene pairs will be evaluated. | |
| --n_batches | Instead of running all of the region x gene pairs, split into this many batches. | |
| --batch | Run this batch (0-indexed) | |
| -parsimony | Logging amount. 1 means that everything is being logged; 10 (the default) is only high-level messages | |


## run_tensorqtl.py

This script takes PLINK-formatted genotypes and parquet-formatted phenotypes, and runs GWAS using the [tensorqtl](https://github.com/broadinstitute/tensorqtl/tree/master/tensorqtl) library.

| Argument | Description | Example Usage |
| --- | --- | --- |
| -plink_genotype | Specify the plink genotype prefix. The script will add the normal `.bed`, `.bim`, and `.fam` file extensions. | ~/genotype/plink_geno_files |
| --parquet_genotype_metadata | If a metadata file for this genotype data exists, this can be merged in and result in more complete GWAS output files | |
| -parquet_phenotype | Specify the path to a parquet phenotype file | | 
| --pval_filter | Only write the results with pvalue less than or equal to this number. Defaults to 1. | | 
| --maf_filter | Only perform GWAS on SNPs with $MAF \geq$ this number. Defaults to 0.05 | |
| -parsimony | Logging amount. 1 means that everything is being logged; 10 (the default) is only high-level messages | |


## validate_model.py

Given the text output of the [model_training.py](#model_trainingpy) script and a held-out data file, this computes prediction performance measures on the heldout data. 


| Argument | Description | Example Usage |
| --- | --- | --- |
| -features | Either a single parquet genotype file, or a parquet genotype file that is configurable with a `chr` key. | ~/genotype/chr-{chr}_variants.parquet |
| -feautes_annotation | Either a single parquet variant metadata file, or a parquet metadata file that is configurable with a `chr` key. | ~/genotype/chr-{chr}_variants_metadata.parquet |
| -data | Phenotype data in parquet format | | 
| -out_fp | Filename to print validation results to. | |
| --model_weights | Text output from the model training script.  | | 
| --sub_batches | Instead of running all of the phenotypes, split into this many batches. | |
| --sub_batch | Run this batch (0-indexed) | |
| -parsimony | Logging amount. 1 means that everything is being logged; 10 (the default) is only high-level messages | |


