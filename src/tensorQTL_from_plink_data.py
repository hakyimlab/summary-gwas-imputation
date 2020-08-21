import sys
import argparse
import logging
import pandas
from genomic_tools_lib.external_tools.tensorqtl import genotypeio, trans
from genomic_tools_lib.external_tools import tensorqtl


# Following the tensorqtl example as closely as possible:
# https://github.com/broadinstitute/tensorqtl/blob/master/example/tensorqtl_examples.ipynb

def run(args):

    # Load phenotypes
    logging.info("Loading phenotypes")
    phenotype_df, phenotype_position_df = tensorqtl.read_phenotype_bed(args.plink_phenotype)
    phenotype_df = phenotype_df.loc[args.idp]

    print(phenotype_df.shape)
    print(phenotype_df.head())

    # Load genotypes
    logging.info("Loading genotypes")
    pr = genotypeio.PlinkReader(args.plink_geno_prefix)
    genotype_df = pr.load_genotypes()

    # Make a fake covariates dataframe
    # TensorQTL requires a covariates dataframe, but an empty one works
    covariates_df = pandas.DataFrame(index=phenotype_df.columns)

    # Run the GWAS function
    logging.info("Running GWAS")
    trans_df = trans.map_trans(genotype_df,
                               phenotype_df,
                               covariates_df,
                               batch_size=2000,
                               pval_threshold=1.0, #return ALL gwas results
                               maf_threshold=0.05)

    # Write output
    logging.info("Writing output")
    trans_df.to_csv(args.out_fp, sep="\t")

    logging.info("Finished")

if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('-plink_phenotype')
    parser.add_argument('-plink_geno_prefix')
    parser.add_argument('-out_fp')
    parser.add_argument('-idp', action='append')

    args = parser.parse_args()

    logging.basicConfig(format='%(asctime)s %(levelname)-8s %(message)s',
                        datefmt='%Y-%m-%d %H:%M:%S',
                        stream=sys.stdout,
                        level=logging.INFO)
    run(args)