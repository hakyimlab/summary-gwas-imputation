import os
import logging
import gzip
from timeit import default_timer as timer
from pyarrow import parquet as pq
import numpy
import pandas
from scipy import stats

from genomic_tools_lib import Logging, Utilities
from genomic_tools_lib.file_formats import Parquet
from genomic_tools_lib.summary_imputation import (Utilities as
                                                  SummaryImputationUtilities,
                                                  SummaryInputation)


GWAS_COLS = ['variant_id', 'panel_variant_id', 'chromosome', 'position',
            'effect_allele', 'non_effect_allele', 'current_build', 'frequency',
            'sample_size', 'zscore', 'pvalue', 'effect_size', 'standard_error',
            'imputation_status', 'n_cases']

# rs554008981	chr1_13550_G_A_b38	chr1	13550	A	G	hg38	0.017316017316017316	337199	-0.9334188062121872	0.35060377415469157	NA	NA	imputed	NA
# rs201055865	chr1_14671_G_C_b38	chr1	14671	C	G	hg38	0.012987012987012988	337199	0.1126679378045526	0.9102938212801959	NA	NA	imputed	NA


#METADATA:
# chromosome: int64
# position: int64
# id: string
# allele_0: string
# allele_1: string
# allele_1_frequency: double
# rsid: string

# FROM MatrixEQTL:
# snps	gene	statistic	pvalue	FDR	beta
# chr10_17966561_C_G	IDP_25029	-12.5922081614529	4.25497460326811e-36	2.13088549409115e-28	-0.185387796514418
# chr10_18221908_T_C	IDP_25029	-12.5579859735497	6.51872918504563e-36	2.13088549409115e-28	-0.180662979423034

def load_matrixEQTL_data(fp):
    """

    :param fp:
    :return:
    """
    df = pandas.read_csv(fp, sep="\t", compression='gzip')
    df.index = df['snps']
    df['standard_error'] = df['beta'] / df['statistic']
    df['zscore'] = df['statistic'] #TODO
    map_dd={'snps': 'panel_variant_id',
            'beta': 'effect_size'}
    df.rename(mapper=map_dd, axis=1, inplace=True)
    return df



def gwas_file_handler(idp, out_dir):
    """

    :param idp:
    :param out_dir:
    :return:
    """
    fname = '{}-gwas-results.txt.gz'.format(idp)
    fp = os.path.join(out_dir, fname)
    if os.path.isfile(fp):
        return gzip.open(fp, 'a')
    else:
        f = gzip.open(fp, 'w')
        header = "\t".join(GWAS_COLS)
        f.write(header.encode())
        return f

def load_parquet_metadata(fp, variants):
    metad_df = pq.read_table(fp).to_pandas()
    metad_df = metad_df.loc[(metad_df['id'].isin(variants))]
    metad_map_dd = {'id': 'variant_id',
                    'allele_0': 'non_effect_allele',
                    'allele_1': 'effect_allele',
                    'allele_1_frequency': 'frequency'}
    metad_df.rename(mapper=metad_map_dd, axis=1, inplace=True)
    return metad_df

def run(args):
    if os.path.exists(args.out):
        logging.error("Output exists. Nope.")
        return
    else:
        os.mkdir(args.out)
    start = timer()
    logging.info("Beginning GWAS conversion")
    for chr in range(1,23):
        matrixeqtl_fp = args.matrixeqtl.format(chr)
        meqtl_df = load_matrixEQTL_data(matrixeqtl_fp)
        logging.log(9, "Loaded matrixEQTL data for chr{}".format(chr))
        variants = set(meqtl_df.index.unique())
        metadata_fp = args.metadata.format(chr)
        metad_df = load_parquet_metadata(metadata_fp, variants)
        logging.log(9, "Loaded metadata for chr{}".format(chr))
        idps = meqtl_df['gene'].unique()
        for i in idps:
            f_i = gwas_file_handler(i, args.out)
            df_i = meqtl_df.join(metad_df, how = 'left')
            print("SHOULD HAVE")
            print(GWAS_COLS)
            print("CURRENTLY HAVE")
            print(list(df_i[GWAS_COLS].columns))

            f_i.close()
            break

        break



if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser("Convert MatrixEQTL results to individual"
                                     "trait GWAS")
    parser.add_argument("-by_region_file",
                        help="If provided, run imputation grouping by regions. "
                             "Much faster.")
    parser.add_argument("-parquet_genotype", help="Parquet Genotype file")
    parser.add_argument("-parquet_genotype_metadata",
                        help="Parquet Genotype variant metadata file")
    parser.add_argument("-gwas_file",
                        help="GWAS file. For the moment, uniform-formatted "
                             "hg38-based files are accepted.")
    parser.add_argument("-window",
                        help="How far to extend in each direction when "
                             "searching for variants", type=int, default=0)
    parser.add_argument("-chromosome", type=int,
                        help="Work only with one chromosome")
    parser.add_argument("-output", help="Where to save stuff")
    parser.add_argument("-cutoff", type=float, default=0,
                        help="naive cutoff when performing SVD")
    parser.add_argument("-regularization", type=float, default=0,
                        help="Ridge-like regularization for matrix inversion")
    parser.add_argument("-frequency_filter", type=float,
                        help="Skip variants with frequency (below f) or "
                             "above (1-f)")
    parser.add_argument("-sub_batches", help="Split the data into subsets",
                        type=int)
    parser.add_argument("-sub_batch", help="only do this subset", type=int)
    parser.add_argument("-containing", help="only do this subset", type=int)
    parser.add_argument("--keep_palindromic_imputation",
                        help="Report imputed values istead of "
                             "(original+flipped) values for palindromic snps",
                        action="store_true")
    parser.add_argument("--use_palindromic_snps",
                        help="Use palindromic variants when imputing",
                        action="store_true")
    parser.add_argument("--standardise_dosages",
                        help="Standardise dosages before computing "
                             "(i.e. use correlation matrix)",
                        action="store_true")
    parser.add_argument("--cache_variants",
                        help="Save variants in memory instead of loading every "
                             "time, when running by variant",
                        action="store_true")
    parser.add_argument("-parsimony",
                        help="Log verbosity level. 1 is everything being "
                             "logged. 10 is only high level messages, above "
                             "10 will hardly log anything",
                        default="10")

    args = parser.parse_args()

    Logging.configure_logging(int(args.parsimony), with_date=True)

    run(args)
