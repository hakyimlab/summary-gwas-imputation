__author__ = "alvaro barbeira"

import gzip
import logging
import re

import numpy
import pandas

from .. import Utilities
from ..Exceptions import ReportableException
from .. import DataSink

def parse_spec(spec, order=None):
    if type(spec) == list:
        order = [x[1] for x in spec] if not order else order
        spec = {x[0]:x[1] for x in spec}
    else:
        raise RuntimeError("Unsupported datframe spec")
    return  spec, order

def load_list(path):
    entries=[]
    for i,line in Utilities.iterate_file(path):
        entries.append(line.strip())
    return entries

def load_column(path, column, unique_entries=True, white_list=None):
    r = set() if unique_entries else list()
    l = (lambda x: r.add(x)) if unique_entries else (lambda x: r.append(x))
    index = None
    for i,line in Utilities.iterate_file(path):
        comps = line.strip().split()

        if i==0:
            index = comps.index(column)
            continue
        v = comps[index]

        if white_list and not v in white_list:
            continue

        l(v)
    return r

def load_dataframe(path, spec=None, order=None, force_special_handling=False, skip_until_header=None, keys=None, key_column_name=None, separator=None, handle_empty_columns=False, additional_filter=None, columns=None):
    #TODO: think of this bash-python command line kink
    if separator == "ANY_WHITESPACE":
        separator = "\s+"

    if force_special_handling or skip_until_header or handle_empty_columns or (keys and key_column_name) or additional_filter:
        d = dataframe_from_text_data_source(path, keys=keys, key_column_name=key_column_name,
                skip_until_header=skip_until_header, separator=separator, handle_empty_columns=handle_empty_columns, additional_filter=additional_filter,
                columns=columns)
    else:
        if separator is None: separator = "\t"
        d = pandas.read_table(path, sep=separator, usecols=columns)

    if spec:
        spec, order = parse_spec(spec, order)
        d = d.rename(columns=spec)
        for c in order:
            if not c in d:
                d= d.assign(**{c:numpy.nan})
        d = d[order]

    return d

#Very much like the previous one but faster, less flexible
def load_dataframe_2(path, keys, key_column_name, spec=None, order=None, to_numeric=None):
    index_column=None
    d = []
    for i, line in Utilities.iterate_file(path):
        if i==0:
            header = line.strip().split()
            if key_column_name and keys:
                index_column = header.index(key_column_name)
            continue
        comps = tuple(line.strip().split())

        if index_column:
            key = comps[index_column]
            if not key in keys:
                continue

        d.append(comps)

    d = Utilities.to_dataframe(data=d, columns=header, to_numeric=to_numeric)

    if spec:
        spec, order = parse_spec(spec, order)
        d = d.rename(columns=spec)
        for c in order:
            if not c in d:
                d= d.assign(**{c:numpy.nan})
        d = d[order]

    return d

def dataframe_from_text_data_source(path, keys=None, key_column_name=None, skip_until_header=None, separator=None, handle_empty_columns=False, sanitize=False, additional_filter=None, columns=None):
    if columns:
        columns = {x for x in columns}

    s = {}
    gz_ = ".gz" in path
    o = gzip.open if gz_ else open
    with o(path) as file:
        header = None
        if skip_until_header:
            _s = skip_until_header if not gz_ else skip_until_header.encode()
            for line in file:
                if _s in line:
                    header = _s
                    c = line.split(_s)
                    if len(c) > 1: header += c[1]
                    break

            if header is None: raise ReportableException("Did not find specified header")
        else:
            header = file.readline()

        if gz_: header = header.decode()

        header_comps = header.strip().split(separator)
        if columns:
            s = {c: [] for c in header_comps if c in columns}
        else:
            s = {c: [] for c in header_comps}

        index = -1
        if key_column_name:
            if not key_column_name in header_comps: raise ReportableException("Did not find key colum name")
            index = header_comps.index(key_column_name)

        header_count = {k:header_comps.count(k) for k in header_comps}
        if len(header_count) < len(header_comps):
            duplicated = [k for k,v in header_count.items() if v>1]
            logging.info("The input GWAS has duplicated columns: %s, will only use the first one in each case", str(duplicated))

        if handle_empty_columns:
            split_r = re.compile(separator) if separator is not None else re.compile("\s")

        for i,line in enumerate(file):
            if gz_: line = line.decode()
            if handle_empty_columns:
                line = line.replace("\n", "")
                comps = split_r.split(line)
            else:
                comps = line.strip().split(separator)

            if sanitize:
                comps = [sanitize_component(x) for x in comps]

            #Yeah, there are those kinds of files
            if not len(comps) == len(header_comps):
                logging.log(8, "Found line with less components than headers, line %i", i)
                continue

            if keys and not comps[index] in keys:
                continue

            if additional_filter and additional_filter(comps):
                continue

            # Load only the first column if in presence of duplicated columns. Yuck!
            sentinel=set()
            for i,c in enumerate(comps):
                comp = header_comps[i]
                if columns and not comp in columns: continue
                if comp in sentinel: continue
                sentinel.add(comp)
                s[comp].append(c)


        for c in header_comps:
            if columns and not c in columns: continue
            s[c] = numpy.array(pandas.to_numeric(s[c], errors='ignore'))

    return pandas.DataFrame(s)

non_en_number = re.compile("^[-\+]?[0-9]*,{1}[0-9]+([eE]{1}[-\+]?[0-9]+)?$")
def sanitize_component(c):
    if non_en_number.match(c): c = c.replace(",",".")
    if c == "": c = None
    elif c == "NA": c = None
    elif c == ".": c = None
    elif c == "\\N": c = None
    elif c == "": c= None
    elif c == "-nan": c = None

    return c

def sanitize_components(comps):
    return [sanitize_component(x) for x in comps]

def to_numeric(d, column):
    if column in d:
        a = [sanitize_component(x) for x in d[column]]
        d[column] = numpy.array(a, dtype=numpy.float64)

#TODO: keep the file object open. Look into pandas code.
class TextDataFrameSink(DataSink.DataFrameSink):
    def __init__(self, path, write_header=True):
        self.path = path
        self.write_header = write_header
        self._wrote_header = False
        self.compression = Utilities._compression(path)

    def sink(self, d):
        header = self.write_header and not self._wrote_header
        mode = "w" if header else "a"
        d.to_csv(self.path, mode=mode, header=header, compression=self.compression, sep="\t", index=False)
        if not self._wrote_header and self.write_header:
            self._wrote_header = True

    def __enter__(self):
        self.initialize()
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.finalize()

    def initialize(self):
        pass

    def finalize(self):
        pass

class TextDataSink(DataSink.DataSink):
    def __init__(self, path, header):
        self.path = path
        self.header = header
        self.file = None

    def initialize(self):
        self.file = gzip.open(self.path, "w")
        self.sink(self.header)

    def finalize(self):
        self.file.close()

    def sink(self, d):
        for _d in d:
            l = "{}\n".format("\t".join(_d)).encode()
            self.file.write(l)

    def __enter__(self):
        self.initialize()
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.finalize()