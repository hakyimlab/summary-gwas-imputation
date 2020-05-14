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
                        "/gpfs/data/im-lab/nas40t2/owen/lab-tools/paper_themes.R"]

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
    def __init__(self, out_dir, in_dir, pheno_name_map):
        """

        :param out_dir: string. Output directory.
        :param in_dir: string. Input directory
        :param pheno_name_map: pandas DataFrame. Columns "pheno_name" and
                                "pheno_num"
        :param load_pheno_f: function. input: pheno_name, pheno_num.
                                function returns a list of p-values for the
                                indicated phenotype
        """
        self.out_dir = out_dir
        self.in_dir = in_dir
        self.pheno_name_map = pheno_name_map
        self.out_name_pattern = os.path.join(self.out_dir,
                                             "metabolite_{}_qq.png")

    def _load_single_pheno(self, pheno_num):
        fp = "chr-{chr}/MatrixEQTL_{num}_chr{chr}.txt"
        fp = os.path.join(self.in_dir, fp)
        all_pvalues = []
        for chr in range(1,23):
            fname_chr = fp.format(num=pheno_num, chr=chr)
            try:
                df = pandas.read_csv(fname_chr, sep = "\t",
                                      usecols=['pvalue'])
                pvalues = list(df['pvalue'])
                all_pvalues.extend(pvalues)
            except FileNotFoundError:
                logging.warning("SS for pheno num {} and chr {} out of place.".format(pheno_num, chr))
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




def load_pheno_map(fp, n_batches=None, batch=None):

    df = pandas.read_csv(fp, sep = "\t")
    ll = list(df.index)
    if (n_batches is not None) and (batch is not None):
        ll = numpy.array_split(ll, n_batches)[batch]

    df = df.loc[ll]
    return df

def run(args):
    Utilities.maybe_create_folder(args.out_dir)
    r_context = RContext()
    name_map = load_pheno_map(args.pheno_map, args.n_batches, args.batch)
    logging.log(9, "Loaded {} phenos".format(len(name_map)))
    file_handler = FileIO(args.out_dir, args.in_dir, name_map)
    for pvals, fp, pheno_name, pheno_num in file_handler.load_all():
        r_context.make_plot(pvals, fp, pheno_name, pheno_num)
    logging.info("Finish Key")


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser()

    parser.add_argument('-pheno_map', help="tsv with columns pheno_name, pheno_num")
    parser.add_argument('-out_dir')
    parser.add_argument('-in_dir')
    parser.add_argument('--n_batches', type=int,
                        help="Break up pheno_map into n_batches chunks")
    parser.add_argument('--batch', type=int, help="Do this chunk. zero-indexed")
    parser.add_argument('--parsimony', type=int, default=9)

    args = parser.parse_args()
    Logging.configure_logging(args.parsimony, target=sys.stdout, with_date=True)
    run(args)

