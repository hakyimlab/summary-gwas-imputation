__author__ = "alvaro barbeira"

import logging
import os
import sqlite3
import pandas
import numpy
import gzip


from genomic_tools_lib import Logging, Utilities
from genomic_tools_lib.data_management import TextFileTools
from genomic_tools_lib.miscellaneous import matrices, Genomics
from genomic_tools_lib.file_formats import Parquet

from covariance_for_model import n_, get_file_map

def run(args):
    if os.path.exists(args.output):
        logging.info("Output already exists, either delete it or move it")
        return

    logging.info("Loading group")
    groups = pandas.read_table(args.group)
    groups = groups.assign(chromosome = groups.gtex_intron_id.str.split(":").str.get(0))
    groups = groups.assign(position=groups.gtex_intron_id.str.split(":").str.get(1))
    groups = Genomics.sort(groups)

    logging.info("Getting parquet genotypes")
    file_map = get_file_map(args)

    logging.info("Getting genes")
    with sqlite3.connect(args.model_db_group_key) as connection:
        # Pay heed to the order. This avoids arbitrariness in sqlite3 loading of results.
        extra = pandas.read_sql("SELECT * FROM EXTRA order by gene", connection)
        extra = extra[extra["n.snps.in.model"] > 0]

    individuals = TextFileTools.load_list(args.individuals) if args.individuals else None

    logging.info("Processing")
    Utilities.ensure_requisite_folders(args.output)

    genes_ = groups[["chromosome", "position", "gene_id"]].drop_duplicates()
    with gzip.open(args.output, "w") as f:
        f.write("GENE RSID1 RSID2 VALUE\n".encode())
        with sqlite3.connect(args.model_db_group_key) as db_group_key:
            with sqlite3.connect(args.model_db_group_values) as db_group_values:
                for i,t_ in enumerate(genes_.itertuples()):
                    g_ = t_.gene_id
                    chr_ = t_.chromosome.split("chr")[1]
                    logging.log(8, "Proccessing %i/%i:%s", i+1, len(genes_), g_)

                    if not n_.search(chr_):
                        logging.log(9, "Unsupported chromosome: %s", chr_)
                        continue
                    dosage = file_map[int(chr_)]

                    group = groups[groups.gene_id == g_]
                    wg=[]
                    for value in group.intron_id:
                        wk = pandas.read_sql("select * from weights where gene = '{}';".format(value), db_group_values)
                        if wk.shape[0] == 0:
                            continue
                        wg.append(wk)

                    if len(wg) > 0:
                        wg = pandas.concat(wg)
                        w = pandas.concat([wk, wg])[["varID", "rsid"]].drop_duplicates()
                    else:
                        w = wk[["varID", "rsid"]].drop_duplicates()

                    if w.shape[0] == 0:
                        logging.log(8, "No data, skipping")
                        continue

                    if individuals:
                        d = Parquet._read(dosage, columns=w.varID.values, specific_individuals=individuals)
                        del d["individual"]
                    else:
                        d = Parquet._read(dosage, columns=w.varID.values, skip_individuals=True)

                    var_ids = list(d.keys())
                    if len(var_ids) == 0:
                        if len(w.varID.values) == 1:
                            logging.log(9, "workaround for single missing genotype at %s", g_)
                            d = {w.varID.values[0]:[0,1]}
                        else:
                            logging.log(9, "No genotype available for %s, skipping",g_)
                            next

                    if args.output_rsids:
                        ids = [x for x in pandas.DataFrame({"varID": var_ids}).merge(w[["varID", "rsid"]], on="varID").rsid.values]
                    else:
                        ids = var_ids

                    c = numpy.cov([d[x] for x in var_ids])
                    c = matrices._flatten_matrix_data([(g_, ids, c)])
                    for entry in c:
                        l = "{} {} {} {}\n".format(entry[0], entry[1], entry[2], entry[3])
                        f.write(l.encode())
    logging.info("Finished building covariance.")

if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser("Generate BSLMM runs on study")
    parser.add_argument("-parquet_genotype_folder", help="Parquet Genotype folder")
    parser.add_argument("-parquet_genotype_pattern", help="Pattern to detect parquet genotypes by chromosome")
    parser.add_argument("-model_db_group_key", help="Model file with group keys as genes")
    parser.add_argument("-model_db_group_values", help="Model file with group values as genes")
    parser.add_argument("-group", help="group definitions")
    parser.add_argument("-output", help="Where to save stuff")
    parser.add_argument("--output_rsids", action="store_true")
    parser.add_argument("--individuals")
    parser.add_argument("-parsimony", help="Log verbosity level. 1 is everything being logged. 10 is only high level messages, above 10 will hardly log anything", default = "10")

    args = parser.parse_args()

    Logging.configure_logging(int(args.parsimony))

    run(args)