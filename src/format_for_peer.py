import os
import pandas as pd
import argparse
import logging
import sys
from genomic_tools_lib import Logging


def load_parquet(fp):
    from pyarrow import parquet as pq
    df = pq.read_table(fp).to_pandas()
    return df


def load_text(fp, sep):
    return pd.read_table(fp, sep=sep)


def write_peer_pheno_df(pheno_df, out_fp):
    """
    PEER requires a DataFrame with no headers or row identifiers.
    This function writes such a file.
    :param pheno_df: DataFrame for writing
    :param out_fp:
    :return: filename
    """
    pheno_df = pheno_df.drop('individual', axis=1)
    pheno_df.to_csv(out_fp, index=False, header=False)
    logging.info("PEER-formatted phenotype file written to {}".format(out_fp))


def run(args):
    in_extension = args.input.split(".")[-1]
    if in_extension == "parquet":
        df = load_parquet(args.input)
    elif in_extension == ".csv":
        df = load_text(args.input, ",")
    elif in_extension == ".txt":
        df = load_text(args.input, "\t")
    else:
        raise ValueError("Input file extension not supported")
    logging.info("Loaded {}".format(args.input))
    write_peer_pheno_df(df, args.output)


if __name__ == "__main__":
    parser = argparse.ArgumentParser("Takes in a phenotype file and creates a"
                                     "PEER-formatted text file")
    parser.add_argument("-input", help="Input filepath")
    parser.add_argument("-output", help="Output filepath")

    args = parser.parse_args()
    Logging.configure_logging(10)

    run(args)





