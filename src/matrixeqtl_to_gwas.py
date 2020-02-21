import os
import logging
import gzip
from timeit import default_timer as timer
from pyarrow import parquet as pq
import pandas

from genomic_tools_lib import Logging


GWAS_COLS = ['variant_id', 'panel_variant_id', 'chromosome', 'position',
            'effect_allele', 'non_effect_allele', 'current_build', 'frequency',
            'sample_size', 'zscore', 'pvalue', 'effect_size', 'standard_error',
            'imputation_status', 'n_cases']

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



def split_file_handler(idp,  out_dir, chromosome):
    """

    :param idp:
    :param out_dir:
    :return:
    """
    fname = '{}_chr-{}_gwas-results.txt.gz'.format(idp, chromosome)
    fp = os.path.join(out_dir, fname)
    return _gwas_file_handler(fp)


def together_file_handler(idp,  out_dir, chromosome):
    """

    :param idp:
    :param out_dir:
    :return:
    """
    fname = '{}_gwas-results.txt.gz'.format(idp)
    fp = os.path.join(out_dir, fname)
    return _gwas_file_handler(fp)


def _gwas_file_handler(fp):
    if os.path.isfile(fp):
        return gzip.open(fp, 'a')
    else:
        f = gzip.open(fp, 'w')
        header = "\t".join(GWAS_COLS) + "\n"
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
    metad_df.index = metad_df['variant_id']
    return metad_df

LINE_STR = "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n"

def gwas_writer(f, df):
    for i in range(len(df)):
        s = LINE_STR.format(*df.iloc[i])
        f.write(s.encode())
    f.close()




def run(args):
    if os.path.exists(args.out):
        logging.error("Output exists. Nope.")
        return
    else:
        os.mkdir(args.out)
    if args.split_by_chromosome:
        file_handler = split_file_handler
    else:
        file_handler = together_file_handler
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
            f_i = file_handler(i, args.out, chr)
            try:
                meqtl_i_df = meqtl_df[meqtl_df['gene']==i]
                df_i = meqtl_i_df.join(metad_df, how = 'left')
                df_i = df_i[GWAS_COLS]
                gwas_writer(f_i, df_i)
            except Exception as e:
                f_i.close()
                raise e
            except KeyboardInterrupt:
                f_i.close()
                logging.error("KEYBOARD INTERRUPT")
                raise KeyboardInterrupt
        logging.log(9, "Completed working with chr{}".format(chr))
    end = timer()
    logging.log(9, "Finished in %.2f seconds" % (end - start))




if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument('-out', help="directory for results")
    parser.add_argument('-matrixeqtl', help="Matrixeqtl filename. Must be able"
                                            " to be formatted with chromosome "
                                            "number.")
    parser.add_argument('-metadata',help="Metadata filename. Must be able "
                                            "to be formatted with chromosome "
                                            "number.")
    parser.add_argument('-parsimony', default=9)
    parser.add_argument('-split_by_chromosome', help="Write an individual file"
                                                     "for each pheno/chromosome"
                                                     "pair.")
    args = parser.parse_args()




    Logging.configure_logging(int(args.parsimony), with_date=True)

    run(args)
