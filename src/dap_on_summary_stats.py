#!/usr/bin/env python
__author__ = "alvaro barbeira"

import os
import logging
import traceback
import re
from collections import namedtuple
from subprocess import call
import pandas
import numpy
import shutil
import pyarrow as pa
from pyarrow import parquet as pq
import time


from timeit import default_timer as timer

from genomic_tools_lib import Logging
from genomic_tools_lib import Utilities
from genomic_tools_lib.Exceptions import ReportableException
from genomic_tools_lib.miscellaneous import PandasHelpers
from genomic_tools_lib.file_formats import Parquet
from genomic_tools_lib.external_tools.dap import Utilities as DAPUtilities, RunDAP


class SummaryStats:
    def __init__(self, meta_path, gene_name, region_re):
        self.meta_path = meta_path
        self.gene_files = self._all_files(meta_path)
        self.gene_name = gene_name
        self.region_re = re.compile(region_re)

    def load_summary_stats(self, summary_stats_path, gene_name):
        d = pandas.read_table(summary_stats_path, usecols=["variant_id",
                                                           "zscore"])
        d = d.rename(columns={"zscore":"z"})
        d["region_id"] = [gene_name] * len(d)
        return d[["region_id", "variant_id", "z"]]

    def _all_files(self, fp):
        with open(fp, 'r') as f:
            s = { i.strip() for i in f.readlines()}
        return s

    def _find_region_name(self, fp):
        s = self.region_re.search(fp.split('/')[-1])
        return s.groups(1)[0]


def _intermediate_folder(intermediate, region): return os.path.join(intermediate, region.region_id)
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

def _run_dap(region, features, features_metadata, summary_stats, intermediate_folder, output_folder, options, dap_command):
    logging.log(9, "fetching and prepatring data")
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

    logging.log(9, "running")
    script_path = _script_path(intermediate_folder, region)
    with open(script_path, "w") as script:
        script.write(command)

    _o = os.path.join(_intermediate_folder(intermediate_folder, region), "dap.o")
    _e = os.path.join(_intermediate_folder(intermediate_folder, region), "dap.e")
    with open(_o, "w") as o:
        with open(_e, "w") as e:
            call(["bash", script_path], stderr=e, stdout=o)
    shutil.move(_o, _output(output_folder, region))
    logging.log(9, "executed dap")

def run_dapg(region, features, features_metadata, summary_stats, intermediate_folder, output_folder, options, dap_command, keep_intermediate=False):
    stats = RunDAP._stats(region.region_id)
    try:
        _run_dap(region, features, features_metadata, summary_stats, intermediate_folder, output_folder, options, dap_command)
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

    return stats



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
    return d

def load_summary_stats(fp, gene_name_re=None, gene_name_col=None):
    if gene_name_re:
        gene_name = _find_gene_name(fp, gene_name_re)

    d = _load(fp, gene_name_col)
    if 'gene_id' in d.columns:
        gene_name = d['gene_id'].iloc[0]
    else:
        gene_name = _find_gene_name(fp, gene_name_re)
        d['gene_id'] = [gene_name] * len(d)
    logging.info("Opening summary stats: {}".format(gene_name))
    return d


def run(args):
    start = timer()

    if os.path.exists(args.output_folder):
        logging.info("Output folder exists. Nope.")
        exit(1)

    os.mkdir(args.output_folder)

    if os.path.exists(args.intermediate_folder):
        logging.warn("Intermediate folder exists.")
        time.sleep(5)
    else:
        os.makedirs(args.intermediate_folder)

    logging.info("Opening features annotation")
    if not args.chromosome:
        features_metadata = pq.read_table(args.parquet_genotype_metadata).to_pandas()
    else:
        features_metadata = pq.ParquetFile(args.parquet_genotype_metadata).read_row_group(args.chromosome-1).to_pandas()

    logging.info("Opening features")
    features = pq.ParquetFile(args.parquet_genotype)

    summary_stats = load_summary_stats(args.summary_stats,
                                       gene_name_re = args.name_re,
                                       gene_name_col=args.gene_col)

    summary_stats = summary_stats[summary_stats.variant_id.isin(features_metadata.id)]
    regions = summary_stats[["region_id"]].drop_duplicates()

    if args.sub_batches is not None and args.sub_batch is not None:
        regions = PandasHelpers.sub_batch(regions, args.sub_batches, args.sub_batch)

    stats = []
    for i, region in enumerate(regions.itertuples()):
        logging.log(9 , "Region %i/%i:%s", i, regions.shape[0], region.region_id)
        _stats = run_dapg(region, features, features_metadata, summary_stats, args.intermediate_folder, args.output_folder, args.options, args.dap_command, not args.keep_intermediate_folder)
        stats.append(_stats)

    stats_path = os.path.join(args.output_folder, "stats.txt")
    stats = RunDAP.data_frame_from_stats(stats).fillna("NA")
    Utilities.save_dataframe(stats, stats_path)

    end = timer()
    logging.info("Ran DAP in %s seconds" % (str(end - start)))



if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser("Generate DAP runs on study")
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
    parser.add_argument("--keep_intermediate_folder", help="don't delete the intermediate stuff", action='store_true')
    parser.add_argument("-parsimony", help="Log verbosity level. 1 is everything being logged. 10 is only high level messages, above 10 will hardly log anything", default = "10")

    args = parser.parse_args()

    Logging.configure_logging(int(args.parsimony))

    run(args)
