#!/usr/bin/env python
__author__ = "alvaro barbeira"

import os
import logging
import traceback
import re
import subprocess
import pandas
import numpy
import shutil
from pyarrow import parquet as pq
import time


from timeit import default_timer as timer

from genomic_tools_lib import Logging
from genomic_tools_lib import Utilities
from genomic_tools_lib.Exceptions import ReportableException
from genomic_tools_lib.miscellaneous import PandasHelpers
from genomic_tools_lib.file_formats import Parquet
from genomic_tools_lib.external_tools.dap import Utilities as DAPUtilities, RunDAP

def _intermediate_folder(intermediate, region): return os.path.join(intermediate, region.region_id)
    # r = region.region_id
    # if type(r) == int:
    #     r = 'region-{}'.format(r)
    # return os.path.join(intermediate, r)
def _stats_path(intermediate, region): return os.path.join(_intermediate_folder(intermediate, region), region.region_id +"_stats.txt")
def _cor_path(intermediate, region): return os.path.join(_intermediate_folder(intermediate, region), region.region_id +"_cor.txt")
def _script_path(intermediate, region): return os.path.join(_intermediate_folder(intermediate, region), region.region_id+".sh")
def _output(output, region): return os.path.join(output, region.region_id+".dap.txt")

r_ = re.compile(r"\\\n[\s]+\\\n")


def _render(s):
    while r_.search(s):
        s = r_.sub("\\\n", s) #substitute empty lines on missing values
    return s

def _dap_command(region, intermediate_folder, output_folder, options, dap_command):
    extra =  " ".join(["-{} {}".format(x[0], x[1]) for x in options]) if options else ""
    command = \
"""#!/usr/bin/env bash
    
[ -d {OUTPUT_DIR} ] || mkdir -p {OUTPUT_DIR}
[ -d {INTERMEDIATE_DIR} ] || mkdir -p {INTERMEDIATE_DIR}
    
{dap} \\
-d_z {stats} \\
-d_ld {cor} \\
{extra}
""".format(OUTPUT_DIR=output_folder,
           INTERMEDIATE_DIR=intermediate_folder,
           dap=dap_command,
           stats=_stats_path(intermediate_folder, region),
           cor=_cor_path(intermediate_folder, region),
           extra=extra)

    command = _render(command)

    return command


def _run_dap(region, features, summary_stats, intermediate_folder,
             output_folder, options, dap_command, snp_f, cluster_f):
    logging.log(9, "fetching and preparing data")
    os.makedirs(_intermediate_folder(intermediate_folder, region))
    s = summary_stats[summary_stats.region_id == region.region_id]
    #m = features_metadata[features_metadata.id.isin(s.variant_id)]
    x = Parquet._read(features, [x for x in s.variant_id.values])
    c = numpy.corrcoef([x[k] for k in s.variant_id.values])
    del x

    numpy.savetxt(_cor_path(intermediate_folder, region), c)
    s[["variant_id", "z"]].rename(columns={"variant_id":"snp_name_i", "z":"z_i"}).to_csv(_stats_path(intermediate_folder, region), index=False, sep="\t")

    del s
    del c

    command = _dap_command(region, intermediate_folder, output_folder, options, dap_command)

    script_path = _script_path(intermediate_folder, region)
    with open(script_path, "w") as script:
        script.write(command)

    dapg_parser = ParseDapGStream(script_path, snp_f, cluster_f)
    logging.log(9, "running")
    dapg_parser.run()
    dapg_parser.write(_output(output_folder, region))
    logging.log(9, "Executed dap")


def run_dapg(region, features, summary_stats, intermediate_folder,
             output_folder, options, dap_command, snp_f, cluster_f,
             keep_intermediate=False):
    """

    :param region:
    :param features: ParquetFile (not yet loaded) of genotype
    :param summary_stats: Pandas DataFrame of summary statistics
    :param intermediate_folder: String. Place to put intermediate files.
    :param output_folder: String. Place to put results.
    :param options: Options for dap-g
    :param dap_command: Path for dap-g command
    :param snp_f: Float.
    :param cluster_f: Float.
    :param keep_intermediate: Boolean.
    :return:
    """
    try:
        _run_dap(region, features, summary_stats, intermediate_folder,
                 output_folder, options, dap_command, snp_f, cluster_f)
    except ReportableException as ex:
        status = Utilities.ERROR_REGEXP.sub('_', ex.msg)
        stats = RunDAP._stats(region.region_id, status=status)
        logging.info("Reportable exception running dap: %s", ex.msg)
    except Exception as ex:
        msg = '{0}'.format(type(ex))  # .replace('"', "'")
        status = '"{0}"'.format(msg)
        stats = RunDAP._stats(region.region_id,  status=status)
        logging.info("Exception running dap:\n%s", traceback.format_exc())
    finally:
        if not keep_intermediate:
            folder = _intermediate_folder(intermediate_folder, region)
            if os.path.exists(folder):
                shutil.rmtree(folder)


class ParseDapGStream:
    def __init__(self, script, snp_filter, cluster_filter):
        self.cmd = ['bash', script]
        self.snp_filter = snp_filter
        self.cluster_filter = cluster_filter
        self.snp_dd = {}

    def run(self):
        logging.log(5, "Running command: {}".format(self.cmd))
        with subprocess.Popen(self.cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE) as proc:
            for line in proc.stdout:
                self._parse(line.decode('utf-8'))

    def _parse(self, s):
        if s.startswith("(("):
            ll = s.split()
            c = int(ll[4])
            if float(ll[2]) >= self.snp_filter and c != -1:
                cluster_lst = self.snp_dd.get(c)
                if cluster_lst is not None:
                    cluster_lst.append("\t".join(ll[1:3]))
                else:
                    self.snp_dd[c] = ["\t".join(ll[1:3])]

        elif s.lstrip().startswith("{"):
            ll = s.split()[:3]
            c = int(ll[0].lstrip("{").rstrip("}"))
            c_pip = float(ll[2])
            if c_pip <= self.cluster_filter:
                val = self.snp_dd.pop(c, None)
                del val

    def write(self, fp):
        header_str = ["snp", "pip", "cluster"]
        with open(fp, 'w') as f:
            f.write("\t".join(header_str) + "\n")
            for k, v in self.snp_dd.items():
                for snp_str in v:
                    f.write(snp_str + "\t" + str(k) + "\n")


def _find_gene_name(fp, regex=None):
    fname = fp.split('/')[-1]

    if regex is None:
        return fname.split('.')[0]
    else:
        regex = re.compile(regex)
        s = regex.search(fname)
        return s.groups(1)[0]


def _load(fp, gene_name_col):
    load_cols = ['variant_id', 'zscore', 'region_id']
    if gene_name_col:
        load_cols.append(gene_name_col)
    d = pandas.read_table(fp, usecols=load_cols)
    map_dd = {gene_name_col: 'gene_id',
              'zscore': 'z'}
    d = d.rename(columns=map_dd)
    d['region_id'] = d['region_id'].astype(str)
    return d


def load_summary_stats(fp, gene_name_re=None, gene_name_col=None):
    
    d = _load(fp, gene_name_col)
    if 'gene_id' in d.columns:
        gene_name = d['gene_id'].iloc[0]
    else:
        gene_name = _find_gene_name(fp, gene_name_re)
        d['gene_id'] = [gene_name] * len(d)
    logging.info("Opening summary stats: {}".format(gene_name))
    return d


def _test_out(dir):
    fname = "test_text.txt"
    content = "Hello world! Printing to {}".format(dir)
    with open(os.path.join(dir, fname), 'w') as f:
        f.write(content)
    logging.log(9, "Tested file output.")


def run(args):
    start = timer()

    if os.path.exists(args.output_folder):
        logging.info("Output folder exists. Nope.")
        exit(1)

    os.mkdir(args.output_folder)

    if os.path.exists(args.intermediate_folder):
        logging.warning("Intermediate folder exists.")
        time.sleep(5)
    else:
        os.makedirs(args.intermediate_folder)

    logging.info("Opening features annotation")
    if not args.chromosome:
        features_metadata = pq.read_table(args.parquet_genotype_metadata,
                                          columns=['id']).to_pandas()
    else:
        features_metadata = pq.ParquetFile(
            args.parquet_genotype_metadata).read_row_group(args.chromosome-1,
                                                           columns=['id']).to_pandas()

    logging.info("Opening features")
    features = pq.ParquetFile(args.parquet_genotype)

    summary_stats = load_summary_stats(args.summary_stats,
                                       gene_name_re = args.name_re,
                                       gene_name_col=args.gene_col)
    if args.testing:
        _test_out(args.output_folder)
        exit(0)

    summary_stats = summary_stats[summary_stats.variant_id.isin(features_metadata.id)]
    regions = summary_stats[["region_id"]].drop_duplicates()

    del features_metadata

    if args.sub_batches is not None and args.sub_batch is not None:
        regions = PandasHelpers.sub_batch(regions, args.sub_batches, args.sub_batch)

    for i, region in enumerate(regions.itertuples()):
        logging.log(9 , "Region %i/%i:%s", i, regions.shape[0], region.region_id)
        run_dapg(region, features, summary_stats, args.intermediate_folder,
                 args.output_folder, args.options, args.dap_command,
                 args.snp_pip, args.cluster_pip,
                 not args.keep_intermediate_folder)


    end = timer()
    logging.info("Ran DAP in %s seconds" % (str(end - start)))



if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser("Generate DAP runs on summary statistics")
    parser.add_argument("-dap_command", help="Which gemma command to run")
    parser.add_argument("-frequency_filter", help="If provided, restrict to variants satisfying MAF criteria", type=float)
    parser.add_argument("-options", help="DAP-G command line options", action="append", nargs=2)
    parser.add_argument("-intermediate_folder", help="Folder to use as scratch space")
    parser.add_argument("-parquet_genotype", help="Parquet Genotype file")
    parser.add_argument("-parquet_genotype_metadata", help="Parquet Genotype variant metadata file")
    parser.add_argument("-summary_stats", help="Summary stats of phenotype being tested")
    parser.add_argument("-name_re", help="Find the phenotype name from the "
                                         "summary stats filename.")
    parser.add_argument("-gene_col", help="Column name of gene if stored in the"
                                          "summary stats.")
    parser.add_argument("-output_folder", help="Where will the model output weights be saved")
    parser.add_argument("-sub_batches", help="Split the data into subsets", type=int)
    parser.add_argument("-sub_batch", help="only do this subset", type=int)
    parser.add_argument("-chromosome", help="Split the data into subsets", type=int)
    parser.add_argument("-snp_pip", help="Keep only the SNPs with PIP "
                                            "greater or equal than this.",
                        type=float)
    parser.add_argument("-cluster_pip", help="Keep only the clusters with PIP "
                                            "greater or equal than this.",
                        type=float)
    parser.add_argument("--keep_intermediate_folder", help="don't delete the intermediate stuff", action='store_true')
    parser.add_argument("-parsimony", help="Log verbosity level. 1 is everything being logged. 10 is only high level messages, above 10 will hardly log anything", default = "10")
    parser.add_argument("--testing", help="Specify to load files, make sure everything is in place, then exit.", default=False, action='store_true')

    args = parser.parse_args()

    Logging.configure_logging(int(args.parsimony))

    run(args)
