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
from genomic_tools_lib import Utilities
from genomic_tools_lib.file_formats import Parquet



class RContext:
    def __init__(self):
        self.sources = [os.path.join(os.path.dirname(__file__), 'generate_qq_plots.R'),
                        "/vol/bmd/meliao/software/other_software/lab-tools/paper_themes.R"]

        self._init_R()

    def _init_R(self):
        R = robjects.r
        pandas2ri.activate()
        logging.log(9, "Initialized R instance.")
        source_cmd = "source('{source}')"
        for i in self.sources:
            R(source_cmd.format(source=i))
        logging.log(9, "Loaded libraries and functions into R.")
        self.R = R
        self._globalenv = robjects.globalenv
        self._make_plot = self.R['make_pval_plot']

    def make_plot(self, pvalues, fp, pheno_name, pheno_num):
        pheno_tag = str(pheno_num) + ": " + pheno_name
        pvals_r = robjects.FloatVector(pvalues)
        logging.log(5, "Beginning R call: {}".format(pheno_tag))
        logging.log(4, "Number of p_values: {}".format(len(pvals_r)))
        self._make_plot(pvals_r, fp, pheno_tag)
        logging.log(5, "Plotted {}".format(pheno_tag))

class FileIO:
    def __init__(self, out_dir, in_dir, pheno_name_map_fp, n_batches, batch):
        """

        :param out_dir: string. Output directory.
        :param in_dir: string. Input directory
        :param pheno_name_map_fp: path to tsv. Should be loaded as a
                                pandas DataFrame. Columns "pheno_name" and
                                "pheno_num"
        :param load_pheno_f: function. input: pheno_name, pheno_num.
                                function returns a list of p-values for the
                                indicated phenotype
        """
        self.out_dir = out_dir
        self.in_dir = in_dir
        self.K_FP = "chr-{chr}/tensorqtl-summ-stats_{pheno}_chr{chr}.txt.gz"
        self.RE_FNAME = re.compile("tensorqtl-summ-stats_(.*)_chr(.*).txt.gz")
        self.K_DIR = "chr-{chr}"
        if pheno_name_map_fp is not None:
            self.pheno_name_map = self._load_pheno_map(pheno_name_map_fp,
                                                       n_batches,
                                                       batch)
        else:
            self.pheno_name_map = self._discover_phenos(n_batches, batch)
        logging.info("Working on {} phenos".format(self.pheno_name_map.shape[0]))
        self.out_name_pattern = os.path.join(self.out_dir,
                                             "ukb-image-gwas_{}_qq.png")

    def _load_single_pheno(self, pheno_num):
        fp = os.path.join(self.in_dir, self.K_FP)
        all_pvalues = []
        for chr in range(1,23):
            fname_chr = fp.format(pheno=pheno_num, chr=chr)
            try:
                df = pandas.read_csv(fname_chr, sep = "\t",
                                      usecols=['pvalue'])
                pvalues = list(df['pvalue'])
                all_pvalues.extend(pvalues)
            except FileNotFoundError:
                logging.warning("SS for pheno {} and chr {} out of place.".format(pheno_num, chr))
            except ValueError:
                dd = pandas.read_csv(fname_chr)
                logging.warning("Bad col names: " + fname_chr)
                logging.warning(' '.join(list(dd.columns)))
        return all_pvalues

    def load_all(self):
        for i, f in self.pheno_name_map.iterrows():
            pvals = self._load_single_pheno(f.pheno_num)
            if len(pvals) == 0:
                continue
            out_fp = self.out_name_pattern.format(f.pheno_num)
            yield (pvals, out_fp, f.pheno_name, f.pheno_num)

    @staticmethod
    def _load_pheno_map(fp, n_batches=None, batch=None):
        df = pandas.read_csv(fp, sep="\t")
        ll = list(df.index)
        if (n_batches is not None) and (batch is not None):
            ll = numpy.array_split(ll, n_batches)[batch]

        df = df.loc[ll]
        return df

    def _discover_phenos(self, n_batches=None, batch=None):
        phenos = set()
        for chr in range(1,23):
            try:
                search_dir = os.path.join(self.in_dir, self.K_DIR.format(chr=chr))
                input_files = os.listdir(search_dir)
                for file in input_files:
                    search_obj = self.RE_FNAME.search(file)
                    if search_obj:
                        phenos.add(search_obj.group(1))
            except FileNotFoundError:
                logging.warning("Results for chromosome {} not found".format(chr))
        df = pandas.DataFrame(data=list(phenos), columns=['pheno_name'])
        df['pheno_num'] = df['pheno_name']
        logging.info("Discovered {} phenotypes".format(df.shape[0]))

        if (n_batches is not None) and (batch is not None):
            ll = list(df.index)
            ll = numpy.array_split(ll, n_batches)[batch]
            df = df.loc[ll]
        return df




def run(args):
    Utilities.maybe_create_folder(args.out_dir)
    r_context = RContext()
    file_handler = FileIO(args.out_dir, args.in_dir, args.pheno_map,
                          args.n_batches, args.batch)
    for pvals, fp, pheno_name, pheno_num in file_handler.load_all():
        r_context.make_plot(pvals, fp, pheno_name, pheno_num)
    logging.info("Finish Key")


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser()

    parser.add_argument('--pheno_map', help="tsv with columns pheno_name, "
                                            "pheno_num")
    parser.add_argument('-out_dir')
    parser.add_argument('-in_dir')
    parser.add_argument('--n_batches', type=int,
                        help="Break up pheno_map into n_batches chunks")
    parser.add_argument('--batch', type=int, help="Do this chunk. zero-indexed")
    parser.add_argument('--parsimony', type=int, default=9)

    args = parser.parse_args()
    Logging.configure_logging(args.parsimony, target=sys.stdout, with_date=True)
    run(args)

