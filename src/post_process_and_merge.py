__author__ = "alvaro barbeira"
import os
import logging

import pandas
import yaml
try:
    from yaml import CLoader as Loader, CDumper as Dumper
except ImportError:
    from yaml import Loader, Dumper

from genomic_tools_lib import Utilities, Logging

def get_trait_map(path):
    with open(path) as f:
        f = yaml.safe_load(f)
        name_map = f["ukb_name"] if "ukb_name" in f else f["name"]
    return name_map

def get_gene_map(path):
    gene_annotation = pandas.read_table(path, usecols=["gene_id", "gene_name"])
    gene_id_map = {}
    gene_name_map = {}
    for t in gene_annotation.itertuples():
        gene_id_map[t.gene_id.split(".")[0]] = t.gene_id
        gene_name_map[t.gene_id] = t.gene_name
    return gene_id_map, gene_name_map

def get_header_names(header):
    if header:
        return header
    return None

def fast_enloc_postprocessing(d, gene_id_map, gene_name_map):
    d = d.assign(gene_id=d.signal_name.apply(lambda x: gene_id_map[x.split(':')[0]]))
    d = d.assign(gene_name=d.gene_id.apply(lambda x: gene_name_map[x]))
    d = d.assign(signal=d.signal_name.apply(lambda x: x.split(':')[1]))
    d = pandas.DataFrame.copy(d)

    d = d.groupby("gene_id").agg(
        gene_name=pandas.NamedAgg(column="gene_name", aggfunc=lambda x: x.values[0]),
        n_signals=pandas.NamedAgg(column="signal", aggfunc=lambda x: x.shape[0]),
        n_snps=pandas.NamedAgg(column="n_snps", aggfunc='sum'),
        eqtl_pip=pandas.NamedAgg(column="eqtl_pip", aggfunc='sum'),
        gwas_pip_woe=pandas.NamedAgg(column="gwas_pip_woe", aggfunc='sum'),
        gwas_pip_we=pandas.NamedAgg(column="gwas_pip_we", aggfunc='sum'),
        rcp=pandas.NamedAgg(column="rcp", aggfunc='sum')
    ).reset_index()
    return d

def run(args):
    if os.path.exists(args.output):
        logging.info("Output exists. Nope.")
        return

    Utilities.ensure_requisite_folders(args.output)
    logging.info("Acquiring files")

    logic = Utilities.file_logic_2(args.input_folder, args.input_pattern, args.name_subfield, args.input_filter)

    trait_map = None
    if args.trait_map:
        logging.info("Loading file mapping")
        trait_map = get_trait_map(args.trait_map)

    gene_id_map, gene_name_map = None, None
    if args.gene_annotation:
        logging.info("Loading gene annotation")
        gene_id_map, gene_name_map = get_gene_map(args.gene_annotation)

    logging.info("Processing files")
    r = []

    for f in logic.itertuples():
        logging.info("Processing %s", f.file)
        names = get_header_names(args.header_names)
        if args.separator == ",":
            d = pandas.read_csv(f.path, header='infer' if not names else None, names=names)
        elif args.separator is None:
            d = pandas.read_table(f.path, header='infer' if not names else None, names=get_header_names(args.header_names), sep="\s+")
        else:
            raise RuntimeError("Unsupported separator")

        if args.specific_post_processing == "FAST_ENLOC":
            d = fast_enloc_postprocessing(d, gene_id_map, gene_name_map)
        elif args.specific_post_processing:
            raise RuntimeError("Unsupported postprocessing option")

        d = d.assign(trait = trait_map[f.trait], tissue = f.tissue)
        r.append(d)

    r = pandas.concat(r)
    if args.integerize_columns:
        logging.info("Convrting columns to integer")
        def integerize_(x):
            try:
                return int(x)
            except:
                return "NA"
        for c in args.integerize_columns:
            r[c] = r[c].apply(integerize_)

    logging.info("Saving")
    Utilities.save_dataframe(r, args.output)

    logging.info("Finished processing.")

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser("Munge fastenloc results into a tidier format")
    parser.add_argument("-input_folder", help="Folder where files lie")
    parser.add_argument("-input_pattern", help="Pattern to parse file names")
    parser.add_argument("--input_filter", help="Pattern to select file names")
    parser.add_argument("-name_subfield", help="Specify multiple key-value pairs to specify filename parsing", nargs=2, action="append", default =[])
    parser.add_argument("-gene_annotation")
    parser.add_argument("-trait_map", help="A file listing conversion for a column")
    parser.add_argument("--header_names", help="If files are headerless, specify column names", nargs="+")
    parser.add_argument("--specific_post_processing", help="optional string asking for a fixed processing to be done on each file. "
                                                          "Currently, 'FAST_ENLOC' is supported")
    parser.add_argument("--separator", help="specific separator for columns")
    parser.add_argument("--integerize_columns", help="forc input columns to have integer values", nargs="+", default=[])
    parser.add_argument("-output", help="Where to save")


    parser.add_argument("-parsimony", help="The higher this value, the less logging you get", type=int, default=10)
    args = parser.parse_args()

    Logging.configure_logging(args.parsimony)
    run(args)