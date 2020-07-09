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
from rpy2.robjects.packages import importr
import subprocess
from genomic_tools_lib import Logging
from genomic_tools_lib.file_formats import Parquet
from genomic_tools_lib.file_formats.gwas import Utilities as GWASUtilities


class RContext:
    def __init__(self, pheno_dd):
        self.RScript = os.path.join(os.path.dirname(__file__), 'matrixeqtl.R')
        self._init_R()
        self.p_mat = self._to_rpy2(pheno_dd)

    def _init_R(self):
        R = robjects.r
        # pandas2ri.activate()
        logging.log(7, "Initialized R instance.")
        _susier = importr('susieR')
        _stats = importr('stats')
        self._susie_rss_f = _susier.susie_rss
        self._susie_pip_f = _susier.susie_get_pip
        self._cor_f = _stats.cor
        logging.log(7, "Loaded libraries and functions into R.")

    @staticmethod
    def _rpy2_to_pandas(df):
        logging.log(5, "Beginning conversion R -> Pandas")
        #        with localconverter(robjects.default_converter + pandas2ri.converter):
        df_pd = pandas2ri.ri2py(df)
        logging.log(5, "Done conversion R -> Pandas")
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
        return y, individuals, names

    @staticmethod
    def _pheno_to_rpy2(pheno, individuals):
        p = pheno.drop_duplicates().reindex(individuals)
        return robjects.FloatVector(p.values)

    def compute_pip(self, geno, ss_generator, fileio, region_id):
        logging.log(5, "Creating R object for genotype")
        print("Geno info: N_individuals: {}, N_varaints: {}".format(len(geno['individual']),
                                                          len(geno) - 1))
        g_mat, individuals, names = self._to_rpy2(geno)
        g_cor = self._cor_f(g_mat)
        ll = len(names)
        pip_df = pandas.DataFrame(columns=['variant_id', 'pip', 'phenotype'])
        for ss_zscores in ss_generator():
            # print(pheno.shape)
            print(ss_zscores.shape)
            if fileio.test_present(region_id, ss_zscores.name):
                logging.log(9,"Skipping pheno already present")
                continue
            z_vec = self._pheno_to_rpy2(ss_zscores, names)
            logging.log(8, "Calculating susie PIP for pheno {}".format(ss_zscores.name))
            susie_fit = self._susie_rss_f(z_vec, g_cor)
            susie_pip = self._susie_pip_f(susie_fit)
            # print(susie_pip[:5])
            # print(numpy.asarray(susie_pip)[:5])
            pips = pandas.DataFrame({'variant_id':names,
                                       'pip': susie_pip,
                                       'phenotype': [ss_zscores.name] * ll})
            fileio.write_region(pips, region_id, ss_zscores.name)


class PythonContext:
    def __init__(self, pheno_fp, metadata_pattern, geno_pattern, ss_dir, ss_format, region_fp, whitelist_fp, sample_size=None, n_batches=None, batch=None):
        """
        0. Save filepaths and file patterns.
        1. Load region_df.
        2. Load (or generate) region x pheno whitelist.
        3. Find intersection of geno, pheno individuals
        4. Load phenotypes
        :param pheno_fp:
        :param metadata_pattern:
        :param geno_pattern:
        :param region_fp:
        :param whitelist_fp:
        :param sample_size:
        :param n_batches:
        :param batch:
        """
        self.GWAS_COLS = ['variant_id', 'panel_variant_id', 'chromosome',
                          'position', 'effect_allele', 'non_effect_allele',
                          'current_build', 'frequency', 'sample_size', 'zscore',
                          'pvalue', 'effect_size', 'standard_error',
                          'imputation_status', 'n_cases', 'gene', 'region_id']
        self._region_df = GWASUtilities.load_region_df(region_fp)
        self._pheno = pq.ParquetFile(pheno_fp)
        self.geno_pattern = geno_pattern
        self.metadata_pattern = metadata_pattern
        self.ss_dir = ss_dir
        self.ss_format = ss_format
        logging.log(3, "Set attributes")
        self.whitelist, phenos = self._load_whitelist(whitelist_fp=whitelist_fp,
                                                          n_batches=n_batches,
                                                          batch=batch)
        self.individuals = self._find_individuals_intersection()
        self.pheno = Parquet._read(self._pheno,
                                   columns=phenos,
                                   specific_individuals=self.individuals)
        self.SAMPLE_SIZE = str(len(self.individuals))
        self._region_df = self._region_df.set_index('region_id')

    def _load_whitelist(self, whitelist_fp=None, n_batches=None, batch=None):
        """
        :param whitelist_fp:
        :param n_batches:
        :param batch:
        :return: tuple (regions, phenos)
                regions: dictionary with keys: region IDs and
                                        values: phenotype names
                phenos: set of phenotype names
        """
        available_phenos = self._pheno.metadata.schema.names
        available_phenos.remove('individual')
        whitelist_dd = {}
        if whitelist_fp is not None:
            logging.info("Loading whitelist provided")
            w_df = pandas.read_table(whitelist_fp)
            w_df = w_df[w_df['gene'].isin(available_phenos)]
            w_df = w_df[w_df['region_id'].isin(self._region_df.region_id)]
            for region, w_df_i in w_df.groupby('region_id'):
                whitelist_dd[region] = set(w_df_i.gene.values)
        else:
            logging.info("Loading whitelist from all pheno x region pairs")
            available_phenos_set = set(available_phenos)
            for region in self._region_df.region_id:
                whitelist_dd[region] = available_phenos_set.copy()

        return self._do_batching(whitelist_dd, n_batches=n_batches, batch=batch)

    @staticmethod
    def _do_batching(dd, n_batches=None, batch=None):
        regions = sorted(list(dd.keys()))
        if (n_batches is not None) and (batch is not None):
            logging.log(9, 'Regions before batching: {}'.format(len(regions)))
            regions = numpy.array_split(regions, n_batches)[batch]
            out_dd = {i:dd[i] for i in regions}
        else:
            out_dd = dd.copy()
        out_phenos = set.union(*[i for i in out_dd.values()])

        logging.info("Working on {} regions".format(len(out_dd)))
        return out_dd, out_phenos

    def _find_individuals_intersection(self):
        pheno_ind_set = self._pheno.read(columns = ['individual'])[0].to_pylist()
        geno_pq = pq.ParquetFile(self.geno_pattern.format(chr=1))
        geno_ind_set = geno_pq.read(columns = ['individual'])[0].to_pylist()

        pheno_ind_set = {str(i) for i in pheno_ind_set}
        geno_ind_set = { str(i) for i in geno_ind_set }

        out_set = pheno_ind_set.intersection(geno_ind_set)
        logging.log(9, "Working with sample size {}".format(len(out_set)))
        return [i for i in out_set]

    def load_region_ss(self, region_name, pheno):
        r_ = self._region_df.loc[region_name]
        print(r_)
        fname_ = self.ss_format.format(chr = r_.chromosome,
                                       pheno = pheno)
        logging.log(3, "Reading {}".format(fname_))
        df = pandas.read_csv(os.path.join(self.ss_dir, fname_), sep="\t",
                             usecols=['panel_variant_id', 'effect_size', 'standard_error', 'position'])
        df = df.loc[(df['position'] < r_.stop) & (df['position'] >= r_.start)]
        df = df.set_index('panel_variant_id')
        df[pheno] = df['effect_size'] / df['standard_error']
        return df[pheno]


    def load_region_genotype(self, region_name):
        region_ = self._region_df.loc[region_name]

        # Load chromosome-specific metadata
        m_fp = self.metadata_pattern.format(chr=region_.chromosome)
        metadata = pq.read_table(m_fp,  columns=['position', 'id']).to_pandas()

        # Find snps in region
        snps = metadata.loc[(metadata.position < region_.stop)
                            & (metadata.position >= region_.start), 'id'].values

        # Load genotype
        g_pq = pq.ParquetFile(self.geno_pattern.format(chr=region_.chromosome))
        return Parquet._read(g_pq,
                             columns=snps,
                             specific_individuals=self.individuals)

    def generate_ss(self, region):
        def f(region=region):
            names = list(self.whitelist[region])
            for name in names:
                yield self.load_region_ss(region, name)
        return f

class FileIO:
    def __init__(self, out_dir, dap_out=True):
        self.out_dir = out_dir
        if dap_out:
            self.fname_format = '{pheno}.txt'
            self.df_preprocess = self._preprocess_dap
        else:
            raise NotImplementedError("Only doing dap output right now.")

    def _filepath(self, r, p):
        return os.path.join(self.out_dir, r, self.fname_format.format(pheno=p))
    def test_present(self, region_id, pheno_name):
        return os.path.isfile(self._filepath(region_id, pheno_name))

    def write_region(self, df, region_id, pheno_name):
        region_dir = os.path.join(self.out_dir, region_id)
        if not os.path.isdir(region_dir):
            os.mkdir(region_dir)
        df_ = self.df_preprocess(df)
        fp = os.path.join(region_dir, self.fname_format.format(pheno=pheno_name))
        df_.to_csv(fp, sep="\t", index=False)

    @staticmethod
    def _preprocess_dap(df):
        df = df.sort_values(by='pip', ascending=False)
        df['region_rank'] = list(range(1, df.shape[0] + 1))
        df['region_rank'] = "((" + df.region_rank.astype(str) + "))"
        return df[['region_rank', 'variant_id', 'pip']]


def run(args):
    if os.path.exists(args.out_dir):
        logging.warning("Output dir already exists at {}".format(args.out_dir))
    else:
        os.mkdir(args.out_dir)

    start_time = timer()
    logging.info("Beginning summary stat calculation")
    f_handler = FileIO(args.out_dir)
    p_context = PythonContext(args.pheno,
                              args.metadata_pattern,
                              args.geno_pattern,
                              args.ss_dir,
                              args.ss_format,
                              args.region,
                              args.whitelist,
                              n_batches=args.n_batches,
                              batch=args.batch)
    r_context = RContext(p_context.pheno)
    n_regions = len(p_context.whitelist)
    for i, region in enumerate(p_context.whitelist.keys()):
        n_phenos_ = len(p_context.whitelist[region])
        logging.log(9, "Working on region {} / {} : {} phenos".format(i + 1, n_regions, n_phenos_))
        geno = p_context.load_region_genotype(region)
        r_context.compute_pip(geno,
                              p_context.generate_ss(region),
                              f_handler,
                              region)
    end_time = timer() - start_time
    logging.info("Finished in {:.2f} seconds".format(end_time))


if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser('Finemapping susieR by region')

    parser.add_argument('-out_dir', help="directory for results", required=True)
    parser.add_argument('-region', help="filepath for regions", required=False)
    parser.add_argument('-metadata_pattern',
                        help="Parquet annotation files, formattable with a {chr} argument",
                        required=True)
    parser.add_argument('-ss_dir', help="Root directory with summary statistics", required=True)
    parser.add_argument('-ss_format', help="Formattable with {chr} and {pheno}", required=True)
    parser.add_argument('-geno_pattern',
                        help="Parquet genotype files, formattable with a {chr} argument",
                        required=True)
    parser.add_argument('-pheno', help="Parquet phenotype file")
    parser.add_argument('--out_split_by', help="What should output be split by?"
                                               " Options are 'region' or "
                                               "'pheno'")
    parser.add_argument('--whitelist', help="Table of region_id x phenotype pairs to run")
    parser.add_argument('--n_batches', type=int)
    parser.add_argument('--batch', type=int)
    parser.add_argument('--parsimony', type=int, default=7)
    parser.add_argument('--compress', help="Gzip all resulting files after writing",
                        default=False, action='store_true')

    args = parser.parse_args()

    Logging.configure_logging(args.parsimony, target=sys.stdout, with_date=True)
    run(args)
