import os
import logging
from timeit import default_timer as timer
from pyarrow import parquet as pq
import pandas
import numpy
import re
import sys
import rpy2.robjects as robjects
from rpy2.robjects import pandas2ri
import subprocess
from genomic_tools_lib import Logging
from genomic_tools_lib.file_formats import Parquet


class RContext:
    def __init__(self, pheno_dd):
        self.RScript = os.path.join(os.path.dirname(__file__), 'matrixeqtl.R')
        self._init_R()
        self.p_mat = self._to_rpy2(pheno_dd)

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
        self._summ_stats_f = self._globalenv['run_summ_stats']

    @staticmethod
    def _rpy2_to_pandas(df):
        logging.log(5, "Beginning conversion R -> Pandas")
        #        with localconverter(robjects.default_converter + pandas2ri.converter):
        df_pd = pandas2ri.ri2py(df)
        logging.log(5, "Finished conversion R -> Pandas")
        return df_pd

    @staticmethod
    def _to_rpy2(dd):
        individuals = dd['individual']
        names = list(dd.keys())
        names.remove('individual')
        x = numpy.array([dd[v] for v in names])
        dimnames = robjects.ListVector([(1, robjects.StrVector(individuals)),
                                        (2, robjects.StrVector(names))])
        y = robjects.r["matrix"](robjects.FloatVector(x.flatten()),
                                 nrow=len(individuals), dimnames=dimnames)
        return y

    def summ_stats(self, geno, i):
        logging.log(5, "Creating R object for genotype")
        g_mat = self._to_rpy2(geno)
        logging.log(9, "Calculating summary statistics for region {}".format(i))
        logging.log(5, "Finished creating string vector. Beginning R call")
        df_r = self._summ_stats_f(g_mat, self.p_mat)
        logging.log(5, "Finished with R call")
        return self._rpy2_to_pandas(df_r)


class PythonContext:
    def __init__(self, pheno_fp, annotation_fp, geno_fp, chr,
                 region_fp=None, sample_size = None, max_r = None,
                 n_batches=None, batch=None):
        # args.pheno, args.annotation, args.chr, args.geno,
        #                               args.region
        self.GWAS_COLS = ['variant_id', 'panel_variant_id', 'chromosome',
                          'position', 'effect_allele', 'non_effect_allele',
                          'current_build', 'frequency', 'sample_size', 'zscore',
                          'pvalue', 'effect_size', 'standard_error',
                          'imputation_status', 'n_cases', 'gene', 'region_id']
        self.chr = int(chr)
        self._pheno = pq.ParquetFile(pheno_fp)
        self.geno = pq.ParquetFile(geno_fp)
        self.individuals = self._find_individuals_intersection()
        self.SAMPLE_SIZE = str(len(self.individuals))
        self.pheno = self._load_phenos(self.individuals, n_batches, batch)
        self.regions = self._load_regions(region_fp, self.chr)
        self.region_index = 0
        self.get_region = self._get_next_region
        self.annotations = self._load_annotations(annotation_fp)
        if max_r is None:
            self.MAX_R = len(self.annotations.region_id.unique())
        else:
            self.MAX_R = max_r
        logging.log(9, "Doing {} regions".format(self.MAX_R))

    def _find_individuals_intersection(self):
        pheno_ind_set = self._pheno.read(columns = ['individual'])[0].to_pylist()
        geno_ind_set = self.geno.read(columns = ['individual'])[0].to_pylist()

        pheno_ind_set = { str(i) for i in pheno_ind_set }
        geno_ind_set = { str(i) for i in geno_ind_set }

        out_set = pheno_ind_set.intersection(geno_ind_set)
        logging.log(9, "Working with sample size {}".format(len(out_set)))
        return [i for i in out_set]

    def _load_phenos(self, individuals, n_batches, batch):
        names = self._pheno.metadata.schema.names
        names.remove('individual')
        if (n_batches is not None) and (batch is not None):
            names = numpy.array_split(names, n_batches)[batch]
        logging.log(9, "Loading data for {} phenos".format(len(names)))
        return Parquet._read(self._pheno, columns=names, specific_individuals=individuals)

    def _load_regions(self, fp, chr):
        if fp is None:
            return None
        else:
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
        cols = list(metad_df.columns)
        cols.append('region_id')
        m_df = pandas.DataFrame(columns = cols)
        if self.regions is not None:
            metad_df['region_id'] = [-1] * l
            for indx, region_ in enumerate(self.regions.itertuples()):
                metad_df.loc[(metad_df['position'] < region_.stop) &
                                      (metad_df['position'] >= region_.start), 'region_id'] = indx
        else:
            metad_df['region_id'] = [0] * len(metad_df)
        metad_df.index = metad_df.variant_id
        logging.log(9, "Loaded metadata")
        return metad_df

    def _get_next_region(self):
        """

        :return: variants, i
        """
        i = self.region_index
        if i is None:
            return None, None
        else:
            variants = self.annotations.loc[self.annotations.region_id == i].variant_id.unique()
            variants = list(variants)
            try:
                g = Parquet._read(self.geno, columns=variants,
                              specific_individuals=self.individuals,
                              to_pandas=False)
            except ValueError:
                print(len(variants))
                var_set = set(variants)
                print(len(var_set))
              #  print([x for x in g.metadata.schema.names])
                exit(1)
            if self.region_index + 1 >= self.MAX_R:
                self.region_index = None
            else:
                self.region_index += 1
            pp = [ x for x in g.keys()]
            if len(pp) > 1:
                return g, i
            else:
                logging.warning("Empty region: {}".format(i))
                return self._get_next_region()

    def to_gwas(self, df):
        df.index = df['snps']
        df['standard_error'] = df['beta'] / df['statistic']
        df['zscore'] = df['statistic']  # TODO
        map_dd = {'snps': 'panel_variant_id',
                  'beta': 'effect_size'}
        df.rename(mapper=map_dd, axis=1, inplace=True)
        df = df.join(self.annotations, how='left')
        # df['region_id'] = [region] * len(df)
        # print(df.columns)
        # print(self.GWAS_COLS)
        return df[self.GWAS_COLS]


class FileIO:
    def __init__(self, out_dir, _split=None):
        self.out_dir = out_dir
        if _split == 'pheno':
            self._fname_format = "MatrixEQTL_{g}_chr{chr}.txt"
            self.K_GROUP = 'gene'
            self.writer = self._write_split
            logging.log(9, "Configured to split region files by pheno")
        elif _split == 'region':
            self._fname_format = "MatrixEQTL_chr{chr}_block{g}.txt"
            self.K_GROUP = 'region_id'
            self.writer = self._write_split
        else:
            self._fname_format = "MatrixEQTL_chr{chr}.txt"
            self.writer = (lambda d, c
                           : self._w(d, self._fname_format.format(chr=c)))

    def _write_split(self, df, chr):
        df_g = df.groupby(self.K_GROUP)
        for name_, g in df_g:
            fp = self._fname_format.format(chr=chr, g = name_)
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


def run(args):
    if os.path.exists(args.out_dir):
        logging.warn("Output dir already exists at {}".format(args.out_dir))
    else:
        os.mkdir(args.out_dir)

    start_time = timer()
    logging.info("Beginning summary stat calculation")
    f_handler = FileIO(args.out_dir, args.out_split_by)
    p_context = PythonContext(args.pheno, args.annotation, args.geno, args.chr,
                              args.region, max_r=args.MAX_R,
                              n_batches=args.n_batches, batch=args.batch)
    r_context = RContext(p_context.pheno)
    ss_df = pandas.DataFrame(columns=['snps', 'gene', 'statistic', 'pvalue', 'FDR', 'beta'])
    while p_context.region_index is not None:
        geno, i = p_context.get_region()
        if i is None:
            continue
        summ_stats = r_context.summ_stats(geno, i)
        ss_df = ss_df.append(summ_stats)
    logging.info("Finished with calculating summary statistics. Beginning file writing")
    gwas_results = p_context.to_gwas(ss_df)
    f_handler.writer(gwas_results, args.chr)
    end_time = timer() - start_time
    logging.log(9, "Finished in {:.2f} seconds".format(end_time))


if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser()

    parser.add_argument('-out_dir', help="directory for results", required=True)
    parser.add_argument('-region', help="filepath for regions", required=False)
    parser.add_argument('-chr', help="Chromosome number", type=int,
                        required=True)
    parser.add_argument('-annotation', help="Parquet annotation file",
                        required=True)
    parser.add_argument('-geno', help="Parquet genotype file")
    parser.add_argument('-pheno', help="Parquet phenotype file")
    parser.add_argument('--out_split_by', help="What should output be split by?"
                                               " Options are 'region' or "
                                               "'pheno'")
    parser.add_argument('--n_batches', type=int)
    parser.add_argument('--batch', type=int)
    parser.add_argument('--parsimony', type=int, default=7)
    parser.add_argument('--compress', help="Gzip all resulting files after writing",
                        default=False, action='store_true')
    parser.add_argument('--MAX_R', type=int)

    args = parser.parse_args()

    Logging.configure_logging(args.parsimony, target=sys.stdout, with_date=True)
    run(args)
