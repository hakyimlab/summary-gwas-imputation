import os
import logging
from timeit import default_timer as timer
from pyarrow import parquet as pq
import pandas
import re

from genomic_tools_lib import Logging

#
# K_SAMPLE_SIZE='10648'
#
# def load_matrixEQTL_data(fp):
#     """
#
#     :param fp:
#     :return:
#     """
#     df = pandas.read_csv(fp, sep="\t", compression='gzip')
#     df.index = df['snps']
#     df['standard_error'] = df['beta'] / df['statistic']
#     df['zscore'] = df['statistic'] #TODO
#     map_dd={'snps': 'panel_variant_id',
#             'beta': 'effect_size'}
#     df.rename(mapper=map_dd, axis=1, inplace=True)
#     return df
#
#
#
# def split_file_handler(idp,  out_dir, chromosome):
#     """
#
#     :param idp:
#     :param out_dir:
#     :return:
#     """
#     fname = '{}_chr-{}_gwas-results.txt.gz'.format(idp, chromosome)
#     fp = os.path.join(out_dir, fname)
#     return _gwas_file_handler(fp)
#
#
# def together_file_handler(idp,  out_dir, chromosome):
#     """
#
#     :param idp:
#     :param out_dir:
#     :return:
#     """
#     fname = '{}_gwas-results.txt.gz'.format(idp)
#     fp = os.path.join(out_dir, fname)
#     return _gwas_file_handler(fp)
#
#
# def _gwas_file_handler(fp):
#     if os.path.isfile(fp):
#         return gzip.open(fp, 'a')
#     else:
#         f = gzip.open(fp, 'w')
#         header = "\t".join(GWAS_COLS) + "\n"
#         f.write(header.encode())
#         return f
#
# def load_parquet_metadata(fp, variants):
#     metad_df = pq.read_table(fp).to_pandas()
#     metad_df = metad_df.loc[(metad_df['id'].isin(variants))]
#     metad_map_dd = {'id': 'variant_id',
#                     'allele_0': 'non_effect_allele',
#                     'allele_1': 'effect_allele',
#                     'allele_1_frequency': 'frequency'}
#     metad_df.rename(mapper=metad_map_dd, axis=1, inplace=True)
#     l = len(metad_df)
#     metad_df['current_build'] = ['hg19'] * l
#     metad_df['sample_size'] = [K_SAMPLE_SIZE] * l
#     metad_df['imputation_status'] = ['NA'] * l
#     metad_df['n_cases'] = ['NA'] * l
#     metad_df.index = metad_df['variant_id']
#     return metad_df
#
# LINE_STR = "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n"
#
# def gwas_writer(f, df):
#     for i in range(len(df)):
#         s = LINE_STR.format(*df.iloc[i])
#         f.write(s.encode())
#     f.close()


class DataHandler:
    def __init__(self, metadata_fp):
        self.K_SAMPLE_SIZE = '10648'
        self.GWAS_COLS = ['variant_id', 'panel_variant_id', 'chromosome',
                          'position', 'effect_allele', 'non_effect_allele',
                          'current_build', 'frequency', 'sample_size', 'zscore',
                          'pvalue', 'effect_size', 'standard_error',
                          'imputation_status', 'n_cases']

        self.metadata_fp = metadata_fp
        self.metadata_df = self._load_parquet_metadata()
        self.loader = self._load_matrixEQTL
        self.writer = self._write_plain

    def _load_parquet_metadata(self):
        metad_df = pq.read_table(self.metadata_fp).to_pandas()
        metad_map_dd = {'id': 'variant_id',
                        'allele_0': 'non_effect_allele',
                        'allele_1': 'effect_allele',
                        'allele_1_frequency': 'frequency'}
        metad_df.rename(mapper=metad_map_dd, axis=1, inplace=True)
        l = len(metad_df)
        metad_df['current_build'] = ['hg19'] * l
        metad_df['sample_size'] = [self.K_SAMPLE_SIZE] * l
        metad_df['imputation_status'] = ['NA'] * l
        metad_df['n_cases'] = ['NA'] * l
        metad_df.index = metad_df['variant_id']
        logging.log(9, "Loaded metadata from {}".format(self.metadata_fp))
        return metad_df

    def _load_matrixEQTL(fp):
        df = pandas.read_csv(fp, sep="\t")
        df.index = df['snps']
        df['standard_error'] = df['beta'] / df['statistic']
        df['zscore'] = df['statistic']  # TODO
        map_dd = {'snps': 'panel_variant_id',
                  'beta': 'effect_size'}
        df.rename(mapper=map_dd, axis=1, inplace=True)
        return df

    def _write_plain(self, df, fp):
        df = df[self.GWAS_COLS]
        df.to_csv(fp, sep="\t", index=False)




class FileHandler:
    def __init__(self, in_dir, out_dir, re_pattern=None, re_flags=None):
        self.in_dir = in_dir
        self.out_dir = out_dir
        # self.file_re = re.compile(re_pattern)
        # self.re_flags = {i[0]: i[1] for i in re_flags}
        self.file_name_df = self._load_file_names(re_pattern, re_flags)

    def _load_file_names(self, re_pattern=None, re_flags=None):
        files_ = os.listdir(self.in_dir)
        logging.log(9, "Loaded files from {}".format(self.in_dir))
        logging.log(9, "Writing files to {}".format(self.out_dir))

        dd = {}

        if re_pattern is not None and re_flags is not None:
            flags_dd = {i[0] : i[1] for i in re_flags}
            file_re = re.compile(re_pattern)
            for i in files_:
                match_i = file_re.search(i)
                for k in flags_dd.keys():
                    dd[k].append(match_i.group(flags_dd[k]))

        dd['input_path'] = [os.path.join(self.in_dir, i) for i in files_]
        dd['output_path'] = [os.path.join(self.out_dir, i) for i in files_]

        return pandas.DataFrame(dd, index=files_)



def convert_files(data_handler, file_handler):
    for _, meqtl_data in file_handler.file_name_df.itertuples():
        meqtl_df = data_handler.loader(meqtl_data.input_path)
        meqtl_df = meqtl_df.join(data_handler.metadata_df, how='left')
        data_handler.writer(meqtl_df, meqtl_data.output_path)




def run(args):
    if os.path.exists(args.out_dir):
        logging.error("Output exists. Nope.")
        return
    else:
        os.mkdir(args.out_dir)

    f_handler = FileHandler(args.in_dir, args.out_dir, args.fname_re,
                            args.fname_field)
    d_handler = DataHandler(args.metadata)

    logging.log("Setup is complete. Beginning GWAS conversion.")
    start = timer()
    convert_files(d_handler, f_handler)
    end = timer()
    logging.log(9, "Finished in %.2f seconds" % (end - start))
    #
    # if args.split_by_chromosome:
    #     file_handler = split_file_handler
    # else:
    #     file_handler = together_file_handler
    # start = timer()
    # logging.info("Beginning GWAS conversion")
    # for chr in range(1,23):
    #     matrixeqtl_fp = args.matrixeqtl.format(chr)
    #     meqtl_df = load_matrixEQTL_data(matrixeqtl_fp)
    #     logging.log(9, "Loaded matrixEQTL data for chr{}".format(chr))
    #     variants = set(meqtl_df.index.unique())
    #     metadata_fp = args.metadata.format(chr)
    #     metad_df = load_parquet_metadata(metadata_fp, variants)
    #     logging.log(9, "Loaded metadata for chr{}".format(chr))
    #     idps = meqtl_df['gene'].unique()
    #     for i in idps:
    #         f_i = file_handler(i, args.out, chr)
    #         try:
    #             meqtl_i_df = meqtl_df[meqtl_df['gene']==i]
    #             df_i = meqtl_i_df.join(metad_df, how = 'left')
    #             df_i = df_i[GWAS_COLS]
    #             gwas_writer(f_i, df_i)
    #         except Exception as e:
    #             f_i.close()
    #             raise e
    #         except KeyboardInterrupt:
    #             f_i.close()
    #             logging.error("KEYBOARD INTERRUPT")
    #             raise KeyboardInterrupt
    #     logging.log(9, "Completed working with chr{}".format(chr))
    # end = timer()
    # logging.log(9, "Finished in %.2f seconds" % (end - start))




if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument('-out_dir', help="directory for results", required=True)
    parser.add_argument('-in_dir', help='directory for output', required=True)
    parser.add_argument('-metadata',help="Metadata filename.", required=True)
    parser.add_argument('-fname_re', help='Regexp pattern for parsing chrom, '
                                          'IDP, and the like.')
    parser.add_argument('-fname_field', help="First arg is group number, "
                                             "second is group name",
                        nargs=2, action='append')

    parser.add_argument('-parsimony', default=9)
    args = parser.parse_args()




    Logging.configure_logging(int(args.parsimony), with_date=True)

    run(args)
