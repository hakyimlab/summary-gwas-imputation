__author__ = "alvaro barbeira"
import pandas

def dapg_signals(file, threshold=None, id_whitelist=None):
    w = pandas.read_table(file, usecols=["gene", "variant_id", "pip", "cluster_id"])
    if id_whitelist:
        w = w[w.variant_id.isin(id_whitelist)]
    w = w.sort_values("pip", ascending=False).groupby(["gene", "cluster_id"]).head(1)
    if threshold:
        w = w[w.pip >= threshold]
    w = w.sort_values("gene")
    return w