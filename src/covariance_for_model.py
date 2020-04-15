__author__ = "alvaro barbeira"

import logging
import os
import re
import sqlite3
import pandas
import numpy
import gzip

from pyarrow import parquet as pq
import pyliftover

from genomic_tools_lib import Logging, Utilities
from genomic_tools_lib.data_management import TextFileTools
from genomic_tools_lib.miscellaneous import matrices, PandasHelpers
from genomic_tools_lib.file_formats import Parquet

def get_file_map(args):
    r = re.compile(args.parquet_genotype_pattern)
    files = os.listdir(args.parquet_genotype_folder)
    files = {int(r.search(f).groups()[0]):os.path.join(args.parquet_genotype_folder, f) for f in files if r.search(f)}
    p = {}
    for k,v in files.items():
        logging.log(4, "Parquet file: {}".format(v))
        g = pq.ParquetFile(v)
        p[k] = g
    return p

def _lift(l, chromosome, position):
    # NA is important, instead of None or NaN, so that integer positions are not converted to floats by pandas. Yuck!
    _new_chromosome = "NA"
    _new_position = "NA"
    try:
        p = int(position)
        l_ = l.convert_coordinate(chromosome, p)
        if l_ is None:
            raise ValueError("Chromosome {} and position {}".format(chromosome, position))
        if l_:
            if len(l_) > 1:
                logging.warning("Liftover with more than one candidate: %s", t.variant_id)
            _new_chromosome = l_[0][0]
            _new_position = int(l_[0][1])
    except:
        pass
    return _new_chromosome, _new_position

def liftover(args, d):
    logging.info("Performing liftover")
    l = pyliftover.LiftOver(args.liftover)
    new_position = []
    new_chromosome = []
    for t in d.itertuples():
        _new_chromosome, _new_position = _lift(l, t.chromosome, t.position)

        new_chromosome.append(_new_chromosome)
        new_position.append(_new_position)

    d = d.assign(chromosome=new_chromosome)
    d = d.assign(position=new_position)
    d = d[d.chromosome.astype(str) !="NA"]
    d = d[d.position.astype(str) != "NA"]

    logging.log(9, "%d variants after liftover", d.shape[0])
    return d

def lift_metadata(metadata, lift_fp):
    l = pyliftover.LiftOver(lift_fp)
    new_position = []
    new_chromosome = []
    for t in metadata.itertuples():
        _new_chrom, _new_pos = _lift(l, t.chromosome, t.position)
        new_chromosome.append(_new_chrom)
        new_position.append(_new_pos)
    metadata = metadata.assign(chromosome = new_chromosome)
    metadata = metadata.assign(position = new_position)

    metadata = metadata[metadata.chromosome.astype(str) != "NA"]
    metadata = metadata[metadata.position.astype(str) != "NA"]
    logging.log(9, "%d variants after liftover", metadata.shape[0])
    metadata = metadata.rename({'id': 'geno_id'}, axis=1)
    metadata = canonical_variant_id(metadata, name='model_id')

    print(metadata.head())
    return metadata[['geno_id', 'model_id']]



def canonical_variant_id(d, name='id', allele_0='allele_0', allele_1='allele_1'):
    cols = ['chromosome', 'position', allele_1, allele_0]
    for i in cols:
        if i not in d:
            raise ValueError("To create variant ID need column: {}".format(i))

    d[name] = (d['chromosome'].astype(str) + '_'
                       + d['position'].astype(str) + "_"
                       + d[allele_0] + "_" + d[allele_1])
    return d

def load_gene_d(dosage, model_variants, geno_variants=None, individuals=None):
    if geno_variants is None:
        geno_variants = model_variants
        g_bool = False
    else:
        logging.log(4, "Using specific geno variants: {}".format(geno_variants[:5]))
        g_bool = True
    if individuals:
        dd = Parquet._read(dosage, columns=geno_variants, specific_individuals=individuals)
        del dd["individual"]
    else:
        dd = Parquet._read(dosage, columns=geno_variants, skip_individuals=True)
    logging.log(4, "Loaded {} snps before renaming".format(len(dd)))
    if g_bool:
        for i in range(len(geno_variants)):
            dd[model_variants[i]] = dd.pop(geno_variants[i])

    return dd



n_ = re.compile("^(\d+)$")

def run(args):
    if os.path.exists(args.output):
        logging.info("Output already exists, either delete it or move it")
        return

    if args.parquet_metadata is not None and args.metadata_lift is not None:
        logging.info("Loading metadata for liftover")
        metadata = pq.read_table(args.parquet_metadata).to_pandas()
        metadata['chromosome'] = 'chr' + metadata['chromosome'].astype(str)
        logging.log(9, "%d variants before liftover", metadata.shape[0])
        lifted_variants = lift_metadata(metadata, args.metadata_lift)
        print(lifted_variants.head())
    else:
        lifted_variants = None

    logging.info("Getting parquet genotypes")
    file_map = get_file_map(args)

    logging.info("Getting genes")
    with sqlite3.connect(args.model_db) as connection:
        # Pay heed to the order. This avoids arbitrariness in sqlite3 loading of results.
        extra = pandas.read_sql("SELECT * FROM EXTRA order by gene", connection)
        #extra = extra[extra["n.snps.in.model"] > 0]

    individuals = TextFileTools.load_list(args.individuals) if args.individuals else None

    logging.info("Processing")
    Utilities.ensure_requisite_folders(args.output)

    with gzip.open(args.output, "w") as f:
        f.write("GENE RSID1 RSID2 VALUE\n".encode())
        with sqlite3.connect(args.model_db) as connection:
            for i,t in enumerate(extra.itertuples()):
                g_ = t.gene
                gene_d = {}
                logging.log(9, "Proccessing %i/%i:%s", i+1, extra.shape[0], g_)
                w = pandas.read_sql("select * from weights where gene = '{}';".format(g_), connection)
                logging.log(5, "Weights loaded for gene {}: {}".format(g_, w.shape[0]))
                if lifted_variants is not None:
                    w = w.merge(lifted_variants, left_on='varID', right_on='model_id')
                    logging.log(5, "Weights after lifted_variants merge: {}".format(w.shape[0]))
                print(w.head())
                w['chr'] = [x[0] for x in w['varID'].str.split("_")]
                w['chr'] = w['chr'].str.lstrip('chr')
                w = w.groupby('chr')
                if len(w.groups.keys()) == 0:
                    logging.log(9, "No chromosomes found for %s, skipping", g_)
                    continue
                for chr_, w_chr in w:
                    if not n_.search(chr_):
                        logging.log(9, "Unsupported chromosome: %s", chr_)
                        continue
                    dosage = file_map[int(chr_)]
                    if 'geno_id' in w_chr:
                        gene_d.update(load_gene_d(dosage,
                                                  w_chr.varID.values,
                                                  geno_variants=w_chr.geno_id.values,
                                                  individuals=individuals))
                    else:
                        gene_d.update(load_gene_d(dosage,
                                                  w_chr.varID.values,
                                                  individuals=individuals))

                var_ids = list(gene_d.keys())
                if len(var_ids) == 0:
                    logging.log(9, "No genotype available for %s, skipping",g_)
                    logging.log(4, 'Could not find {}'.format(w_chr.varID.values[:5]))
                    continue

                if args.output_rsids:
                    ids = [x for x in pandas.DataFrame({"varID": var_ids}).merge(w_chr[["varID", "rsid"]], on="varID").rsid.values]
                else:
                    ids = var_ids

                c = numpy.cov([gene_d[x] for x in var_ids])
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
    parser.add_argument("-model_db", help="Where to save stuff")
    parser.add_argument("-output", help="Where to save stuff")
    parser.add_argument("--output_rsids", action="store_true")
    parser.add_argument("--metadata_lift")
    parser.add_argument("--parquet_metadata")
    parser.add_argument("--individuals")
    parser.add_argument("-parsimony", help="Log verbosity level. 1 is everything being logged. 10 is only high level messages, above 10 will hardly log anything", default = "10")

    args = parser.parse_args()

    Logging.configure_logging(int(args.parsimony))

    run(args)
