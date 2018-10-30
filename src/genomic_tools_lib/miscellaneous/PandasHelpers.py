__author__ = "alvaro barbeira"
import numpy

def sub_batch(d, sub_batches, sub_batch):
    batches = numpy.array_split(range(0, d.shape[0]), sub_batches)
    batch = batches[sub_batch]
    d = d.iloc[batch]
    d = d.reset_index(drop=True)
    return d

def split_column(d, split_column_args):
    column = split_column_args[0]
    separator = split_column_args[1]
    targets = split_column_args[2:]
    c = list(zip(*[x.split(separator) for x in d[column]]))
    args = {t:c[i] for i,t in enumerate(targets)}
    d = d.assign(**args)
    return d
