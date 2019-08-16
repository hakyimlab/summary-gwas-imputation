#!/usr/bin/env python
__author__ = "alvaro barbeira"
import logging
import re
import gzip

import pandas

from genomic_tools_lib import Utilities, Logging

def get_gene_whitelist(args):
    g = pandas.read_table(args.gencode)
    return {x for x in g.gene_id.values}

model_re = re.compile("\s+(?P<rank>\d+)\s+(?P<mpp>[-+]?\d*\.\d+e[-+]?\d+)\s+(?P<n>\d+)\s+(?P<posterior_score>[-+]?\d*\.\d+|\d+)\s+\[(?P<variants>.*)\]\n$")
def parse_model_line(name, s):
    return "{}\t{}\t{}\t{}\t{}\n".format(name, s.group("rank"), s.group("n"), s.group("mpp"), s.group("posterior_score"))

def parse_model_line_for_variant(name, s):
    v = s.group("variants").split("] [")
    if v[0] == "NULL":
        return None
    return ["{}\t{}\t{}\n".format(name, s.group("rank"), x) for x in v]

model_expected_size_re = re.compile("Posterior expected model size: (?P<p>[-+]?\d*\.\d+|\d+) \(sd = (?P<pse>[-+]?\d*\.\d+|\d+)\)\n$")
def parse_expected_size(s):
    return s.group("p"), s.group("pse")

lognc_re = re.compile("LogNC = (?P<lognc>[-+]?\d*\.\d+|\d+) \( Log10NC = (?P<log10nc>[-+]?\d*\.\d+|\d+) \)\n$")
def parse_log_10_nc(s):
    return s.group("lognc"), s.group("log10nc")

variant_re = re.compile("\(\((?P<variant_rank>\d+)\)\)\s+(?P<variant_id>[^\s]+)\s+(?P<variant_pip>[-+]?\d*\.\d+e[-+]?\d+)\s+(?P<variant_log10abvf>[-+]?\d*\.\d+|\d+)\s+(?P<cluster_id>\d+)\n$")
def parse_variant_line(s):
    return s.group("variant_rank"), s.group("variant_id"), s.group("variant_pip"), s.group("variant_log10abvf"), s.group("cluster_id")

cluster_re = re.compile("\s+\{(?P<cluster_id>\d+)\}\s+(?P<n>\d+)\s+(?P<cluster_pip>[-+]?\d*\.\d+e[-+]?\d+)\s+(?P<r2>[-+]?\d*\.\d+|\d+)\s+(?P<correlation>.*)\n$")
def parse_cluster_line(s):
    return s.group("cluster_id"), s.group("n"), s.group("cluster_pip"), s.group("r2")

def run(args):
    logging.info("Processing...")
    Utilities.ensure_requisite_folders(args.output_prefix)

    spec = Utilities.file_logic(args.input_folder, args.input_pattern)

    with gzip.open(args.output_prefix+ ".models.txt.gz", mode="w") as models:
        models.write("gene\tmodel\tn\tpp\tps\n".encode())
        with gzip.open(args.output_prefix + ".models_variants.txt.gz", mode="w") as model_variants:
            model_variants.write("gene\tmodel\tvariant\n".encode())
            with gzip.open(args.output_prefix+ ".model_summary.txt.gz", mode="w") as model_summary:
                model_summary.write("gene\tpes\tpes_se\tlog_nc\tlog10_nc\n".encode())
                with gzip.open(args.output_prefix+ ".variants_pip.txt.gz", mode="w") as variant_pip:
                    variant_pip.write("gene\trank\tvariant_id\tpip\tlog10_abf\tcluster_id\n".encode())
                    with gzip.open(args.output_prefix + ".clusters.txt.gz", mode="w") as clusters:
                        clusters.write("gene\tcluster\tn_snps\tpip\taverage_r2\n".encode())
                        with gzip.open(args.output_prefix + ".cluster_correlations.txt.gz", mode="w") as cluster_correlations:
                            cluster_correlations.write("gene\tid1\tid2\tvalue\n".encode())
                            for i,t in enumerate(spec.itertuples()):
                                logging.log(9, "Processing %s", t.name)
                                written = set()
                                with open(t.path) as dap:
                                    p, pse, lognc, log10nc = None, None, None, None
                                    for l in dap:
                                        s = model_re.search(l)
                                        if s:
                                            ml = parse_model_line(t.name, s)
                                            models.write(ml.encode())
                                            vl = parse_model_line_for_variant(t.name, s)
                                            if vl:
                                                for vl_ in vl:
                                                    model_variants.write(vl_.encode())
                                            continue

                                        s = model_expected_size_re.search(l)
                                        if s:
                                            p, pse = parse_expected_size(s)
                                            continue

                                        s = lognc_re.search(l)
                                        if s:
                                            lognc, log10nc = parse_log_10_nc(s)
                                            model_summary.write("{}\t{}\t{}\t{}\t{}\n".format(t.name, p, pse, lognc, log10nc).encode())
                                            continue

                                        s = variant_re.search(l)
                                        if s:
                                            rank, id, pip, log10_abvf, cluster_id = parse_variant_line(s)
                                            variant_pip.write("{}\t{}\t{}\t{}\t{}\t{}\n".format(t.name, rank, id, pip, log10_abvf, cluster_id).encode())
                                            continue

                                        s = cluster_re.search(l)
                                        if s:
                                            id, n, pip, r2 = parse_cluster_line(s)
                                            clusters.write("{}\t{}\t{}\t{}\t{}\n".format(t.name, id, n, pip, r2).encode())

                                            _id1 = int(id)
                                            comps = s.group("correlation").strip().split()

                                            for _id2 in range(1, len(comps)+1):
                                                if (_id1,_id2) in written or (_id2,_id1) in written:
                                                    continue
                                                comp = comps[_id2-1]
                                                cluster_correlations.write("{}\t{}\t{}\t{}\n".format(t.name, _id1, _id2, comp).encode())
                                                written.add((_id1,_id2))
    logging.info("Finished")

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser("Patch gene name column")
    parser.add_argument("-input_folder", help="Where to load file from")
    parser.add_argument("-input_pattern", help="Name pattern")
    parser.add_argument("-output_prefix", help="Where to save")
    parser.add_argument("-parsimony", help="Logging parsimony", type=int, default=10)
    args = parser.parse_args()

    Logging.configure_logging(args.parsimony)
    run(args)