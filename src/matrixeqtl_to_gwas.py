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
    l = len(metad_df)
    metad_df['current_build'] = ['hg19'] * l
    metad_df['sample_size'] = ['10648'] * l
    metad_df['imputation_status'] = ['NA'] * l
    metad_df['n_cases'] = ['NA'] * l
    return metad_df

LINE_STR = "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n"

def gwas_writer(f, df):
    for i in range(len(df)):
        f.write(LINE_STR.format(*df.iloc[i]).encode())
    f.close()




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
            try:
                df_i = meqtl_df.join(metad_df, how = 'left')[GWAS_COLS]
                gwas_writer(f_i, df_i)
            except Exception as e:
                f_i.close()
                raise e
        logging.log(9, "Completed working with chr{}".format(chr))
        break
    end = timer()
    logging.log(9, "Finished in %.2f seconds" % (start - end))




if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument('-out', default='/vol/bmd/meliao/data/imageXcan/test/meqtl')
    parser.add_argument('-matrixeqtl',
                   default='/vol/bmd/meliao/data/imageXcan/matrixEQTL/2020-02-17_chr-{}_matrixeqtl-out.txt.gz')
    parser.add_argument('-metadata',
                   default='/vol/bmd/meliao/data/parquet/parquet/ukb_imp_chr{}_v3.variants_metadata.parquet')
    args = parser.parse_args()




    Logging.configure_logging(int(args.parsimony), with_date=True)

    run(args)
