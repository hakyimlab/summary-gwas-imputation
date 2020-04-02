import os
import logging
from timeit import default_timer as timer
from pyarrow import parquet as pq
import pandas
import re
import sys
import rpy2.robjects as robjects
from rpy2.robjects import pandas2ri
import subprocess
from genomic_tools_lib import Logging


class RContext:
    def __init__(self, pheno_fp, annotation_fp, chr, geno_fp):
        self.RScript = os.path.join(os.path.dirname(__file__), 'matrix_eQTL.R')
        logging.log(5, "Pheno fp: {}".format(pheno_fp))
        logging.log(5, "Geno fp: {}".format(geno_fp))
        logging.log(5, "Annotation fp: {}".format(annotation_fp))
        self._init_R()
        self.init_lst = self._create_init_data(pheno_fp, annotation_fp, chr,
                                               geno_fp)


    def _init_R(self):
        R = robjects.r
        pandas2ri.activate()
        logging.log(9, "Initialized R instance.")
        source_cmd = "source('{source}')"
        cmd = source_cmd.format(source=self.RScript)
        logging.log(6, "Running R command: {}".format(cmd))
        R(cmd)
        logging.log(9, "Loaded libraries and functions into R.")
        self.R = R
        self._globalenv = robjects.globalenv
        self._init_f = self._globalenv['init_data']
        self._summ_stats_f = self._globalenv['run_summ_stats']


    def _create_init_data(self, pheno, annotation, chr, geno):
        logging.log(9, "Initializing data in R")
        init_lst = self._init_f(pheno, chr, geno)
        logging.log(9, "R data has been initialized.")
        return init_lst

    def _rpy2_to_pandas(self, df):
        logging.log(5, "Beginning conversion R -> Pandas")
#        with localconverter(robjects.default_converter + pandas2ri.converter):
        df_pd = robjects.conversion.rpy2py(df)
        logging.log(5, "Finished conversion R -> Pandas")
        return df_pd


    def summ_stats(self, variants, i):
        logging.log(9, "Calculating summary statistics for region {}".format(i))
        logging.log(5, "Creating a string vector of variants")
        v_ = robjects.StrVector(variants)
        logging.log(5, "Finished creating string vector. Beginning R call")
        df_r = self._summ_stats_f(self.init_lst, v_)
        logging.log(5, "Finished with R call")
        logging.log(9, "Done calculating summary statistics for region %i " % i)
        return self._rpy2_to_pandas(df_r)

class PythonContext:
    def __init__(self, region_fp, chr, annotation_fp):
        self.SAMPLE_SIZE = '908'
        self.GWAS_COLS = ['variant_id', 'panel_variant_id', 'chromosome',
                          'position', 'effect_allele', 'non_effect_allele',
                          'current_build', 'frequency', 'sample_size', 'zscore',
                          'pvalue', 'effect_size', 'standard_error',
                          'imputation_status', 'n_cases', 'gene', 'region_id']
        self.chr = int(chr)
        self.regions = self._load_regions(region_fp, self.chr)
        self.region_index = 0
        self.get_region = self._get_next_region
        self.annotations = self._load_annotations(annotation_fp)

    def _load_regions(self, fp, chr):
        region_df = pandas.read_table(fp)
        region_df.columns = [i.strip() for i in region_df.columns]
        region_df['chromosome'] = region_df['chr'].str.lstrip('chr').astype(int)
        region_df = region_df.loc[region_df['chromosome'] == chr]
        logging.log(9, "Loaded {} regions".format(len(region_df)))

        return region_df

    def _load_annotations(self, fp):
        metad_df = pq.read_table(fp).to_pandas()
        metad_map_dd = {'id': 'variant_id',
                        'allele_0': 'non_effect_allele',
                        'allele_1': 'effect_allele',
                        'allele_1_frequency': 'frequency'}
        metad_df.rename(mapper=metad_map_dd, axis=1, inplace=True)
        l = len(metad_df)
        metad_df['current_build'] = ['hg19'] * l
        metad_df['sample_size'] = [self.SAMPLE_SIZE] * l
        metad_df['imputation_status'] = ['NA'] * l
        metad_df['n_cases'] = ['NA'] * l
        metad_df.index = metad_df['variant_id']
        logging.log(9, "Loaded metadata")
        logging.log(5, "Loaded metadata from {}".format(fp))
        return metad_df

    def _get_next_region(self):
        """

        :return: variants, i
        """
        i = self.region_index
        if i is None:
            raise ValueError("get_region() should not be queried now ")
        else:
            start = self.regions.iloc[i].start
            stop = self.regions.iloc[i].stop
            annotations_i = self.annotations.loc[(start <= self.annotations.position) &
                                                    (self.annotations.position <=
                                                 stop)]
            variants = annotations_i.variant_id
            if self.region_index + 1 == len(self.regions):
                self.region_index = None
            else:
                self.region_index += 1
            return variants, i

    def to_gwas(self, df, region):
        df.index = df['snps']
        df['standard_error'] = df['beta'] / df['statistic']
        df['zscore'] = df['statistic']  # TODO
        map_dd = {'snps': 'panel_variant_id',
                  'beta': 'effect_size'}
        df.rename(mapper=map_dd, axis=1, inplace=True)
        df = df.join(self.annotations, how = 'left')
        df['region_id'] = [region] * len(df)
        # print(df.columns)
        # print(self.GWAS_COLS)
        return df[self.GWAS_COLS]


class FileIO:
    def __init__(self, out_dir, split_by_pheno=False):
        self.out_dir = out_dir
        if split_by_pheno:
            self._fname_format = "MatrixEQTL_{pheno}_chr{chr}_block{index}.txt"
            self.K_GROUP = 'gene'
            self.writer = self._write_split_pheno
            logging.log(9, "Configured to split region files by pheno")
        else:
            self._fname_format = "MatrixEQTL_chr{chr}_block{index}.txt"
            self.writer = self._write_region
            logging.log(9, "Configured to write whole region into one file")

    def _write_split_pheno(self, df, chr, indx):
        df_g = df.groupby(self.K_GROUP)
        for pheno, g in df_g:
            fp = self._fname_format.format(pheno=pheno, chr=chr, index=indx)
            self._w(g, fp)
        pass

    def _write_region(self, df, chr, indx):
        fp = self._fname_format.format(chr=chr, index=indx)
        self._w(df, fp)

    def _w(self, df, fp):
        out_fp = os.path.join(self.out_dir, fp)
        logging.log(3, "Writing file to {}".format(out_fp))
        df.to_csv(out_fp, sep="\t", index=False)

    def compressor(self, dir=None):
        if dir is None:
            dir = self.out_dir
        logging.log(9, "Compressing all files in {}".format(dir))
        cmd = ['gzip', '-r', dir]
        subprocess.run(cmd)
        logging.log(9, "Finished with compression")




def convert_files(data_handler, file_handler):
    for meqtl_data in file_handler.file_name_df.itertuples():
        meqtl_df = data_handler.loader(meqtl_data.input_path)
        meqtl_df = meqtl_df.join(data_handler.metadata_df, how='left')
        data_handler.writer(meqtl_df, meqtl_data.output_path)




def run(args):
    if os.path.exists(args.out_dir):
        logging.error("Output exists. Nope.")
        return
    else:
        os.mkdir(args.out_dir)

    start_time = timer()
    logging.log(9, "Beginning summary stat calculation")
    f_handler = FileIO(args.out_dir, args.split_by_pheno)
    r_context = RContext(args.pheno, args.annotation, args.chr, args.geno)
    p_context = PythonContext(args.region, args.chr, args.annotation)
    while p_context.region_index is not None:
        variants, i = p_context.get_region()
        summ_stats = r_context.summ_stats(variants, i)
        gwas_results = p_context.to_gwas(summ_stats, i)
        f_handler.writer(gwas_results, args.chr, i)
    end_time = timer() - start_time
    logging.log(9, "Finished in {:.2f} seconds".format(end_time))




if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser()

    parser.add_argument('-out_dir', help="directory for results", required=True)
    parser.add_argument('-region', help="filepath for regions", required=True)
    parser.add_argument('-chr', help="Chromosome number", type=int,
                        required=True)
    parser.add_argument('-annotation', help="Parquet annotation file",
                        required=True)
    parser.add_argument('-geno', help="Parquet genotype file")
    parser.add_argument('-pheno')
    parser.add_argument('--split_by_pheno', default=False, action='store_true')
    parser.add_argument('--parsimony', type=int, default=7)
    parser.add_argument('--compress', "Gzip all resulting files after writing",
                        default=False, action='store_true')


    args = parser.parse_args()

    Logging.configure_logging(args.parsimony, target=sys.stdout, with_date=True)
    run(args)
