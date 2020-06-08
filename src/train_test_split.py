import logging
from timeit import default_timer as timer
import pandas as pd
import numpy as np
import sys
from genomic_tools_lib import Logging

def load_parquet(fp):
    from pyarrow import parquet as pq
    df = pq.read_table(fp).to_pandas()
    return df


def load_text(fp, sep):
    return pd.read_table(fp, sep=sep)


def find_train_ids(ids, train_frac, n_train, seed):
    ll = len(ids)
    if n_train and ll > n_train:
        n_ = n_train
    else:
        n_ = np.floor(ll * train_frac)
    np.random.seed(seed)
    train_ids = np.random.choice(ids, int(n_), replace=False)
    return set(train_ids)

def write_parquet(df, prefix):
    fp = prefix + ".parquet"
    df.to_parquet(fp)
    logging.info("Wrote to: {}".format(fp))

def write_text(df, prefix):
    fp = prefix + ".txt"
    df.to_csv(fp, sep="\t")
    logging.info("Wrote to: {}".format(fp))



def run(args):
    start = timer()
    in_extension = args.input.split(".")[-1]
    if in_extension == "parquet":
        df = load_parquet(args.input)
    elif in_extension == "csv":
        df = load_text(args.input, ",")
    elif in_extension == "txt":
        df = load_text(args.input, "\t")
    else:
        raise ValueError("Input file extension not supported")
    logging.info("Loaded {}".format(args.input))
    train_ids = find_train_ids(df.individual.values, args.train_frac,
                                   args.n_train, args.seed)
    logging.info("Length of loaded dataset: {}".format(df.shape[0]))
    train_df = df[df['individual'].isin(train_ids)]
    logging.info("Length of train dataset: {}".format(train_df.shape[0]))
    test_df = df[~df['individual'].isin(train_ids)]
    logging.info("Length of test dataset: {}".format(test_df.shape[0]))
    train_pre = args.output_prefix + "_train"
    test_pre = args.output_prefix + "_test"

    write_parquet(train_df, train_pre)
    write_text(train_df, train_pre)

    write_parquet(test_df, test_pre)
    write_text(test_df, test_pre)
    end = timer() - start
    logging.info("Finished in {:.2f} seconds".format(end))




if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser()

    parser.add_argument('-input', help="Input phenotype file", required=True)
    parser.add_argument('-output_prefix')
    parser.add_argument('-seed', type=int, default=10000)
    parser.add_argument('-train_frac', default=0.8, type=float)
    parser.add_argument('-n_train', type=int)
    parser.add_argument('--parsimony', type=int, default=7)

    args = parser.parse_args()

    Logging.configure_logging(args.parsimony, target=sys.stdout, with_date=True)
    run(args)
