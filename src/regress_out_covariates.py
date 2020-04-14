__author__ = "alvaro barbeira"
import logging

import pandas
from scipy import stats
from patsy import dmatrices
import statsmodels.api as sm
import pyarrow as pa
import pyarrow.parquet as pq

from genomic_tools_lib import Logging, Utilities
from genomic_tools_lib.file_formats import Parquet

def inverse_normalize(df):
    """
    Takes a DataFrame, draws a random Normal sample, and quantile normalizes
    the series to the reference Normal sample.
    :return: pandas.DataFrame
    """
    rank_df = df.rank(method='average')
    rank_df = rank_df / (len(rank_df) + 1)
    quantiles = stats.norm.ppf(rank_df)
    return pandas.DataFrame(quantiles, index=df.index, columns=df.columns)

def run(args):
    logging.info("Starting")
    Utilities.ensure_requisite_folders(args.output)

    logging.info("Read covariate")
    covariate = pq.read_table(args.covariate).to_pandas()
    logging.info("Read data")
    data = pq.read_table(args.data).to_pandas()

    if args.inverse_normalize_data:
        data = inverse_normalize(data)
    logging.info("Processing")
    covariate_names = covariate.columns.values[1:]
    results = {"individual":data.individual.values}
    variables = [x for x in data.columns.values[1:]]
    for i,column in enumerate(variables):
        logging.log(9, "%i/%i:%s", i, len(variables), column)
        d = data[["individual", column]].rename(columns={column:"y"}).merge(covariate, on="individual", how="inner").drop("individual", axis=1)
        y, X = dmatrices("y ~ {}".format(" + ".join(covariate_names)), data=d, return_type="dataframe")
        model = sm.OLS(y, X)
        result = model.fit()
        results[column] = result.resid
    results = pandas.DataFrame(results)[["individual"]+variables]
    Parquet.save_variable(args.output, results)
    logging.info("Finished")

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser("Get data residuals of the linear fit of specific covariates")
    parser.add_argument("-covariate")
    parser.add_argument("-data")
    parser.add_argument("--inverse_normalize_data", default=False,
                        action='store_true')
    parser.add_argument("-output")
    parser.add_argument("-parsimony", type=int, default=10)
    args = parser.parse_args()
    Logging.configure_logging(args.parsimony)
    run(args)