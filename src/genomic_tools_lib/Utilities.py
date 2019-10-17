__author__ = "alvaro barbeira"

import os
import re
import gzip
import logging
import pandas
import csv

from .Exceptions import ReportableException

def folder_contents(folder, pattern=None):
    regexp = re.compile(pattern) if pattern else None
    p = os .listdir(folder)
    if regexp: p = [x for x in p if regexp.search(x)]
    p = [os.path.join(folder,x) for x in p]
    return p

def ensure_requisite_folders(path):
    folder = os.path.split(path)[0]
    if len(folder) and not os.path.exists(folder):
        os.makedirs(folder)

def maybe_create_folder(path):
    if not os.path.exists(path):
        os.makedirs(path)

########################################################################################################################

def file_logic(folder, pattern):
    r = re.compile(pattern)
    f = sorted([x for x in os.listdir(folder) if r.search(x)])
    p_ = [r.search(x).group(1) for x in f]
    p = [os.path.join(folder,x) for x in f]
    return pandas.DataFrame({"name":p_, "path":p, "file":f})

def file_logic_2(folder, pattern, sub_fields, filter=None):
    r = re.compile(pattern)

    f = os.listdir(folder)
    if filter is not None:
        filter = re.compile(filter)
        f = [x for x in f if filter.search(x)]

    f = sorted([x for x in f if r.search(x)])
    r, subfield_names, subfield_positions = name_parse_prepare(pattern, sub_fields)
    values=[]
    for x in f:
        values.append((x, os.path.join(folder, x))+name_parse(x, r, subfield_positions))
    columns = ["file", "path"] + subfield_names
    values = pandas.DataFrame(values, columns=columns)
    return values

########################################################################################################################

def name_parse_prepare(name_subfield_regexp, name_subfield):
    if name_subfield and name_subfield_regexp:
        r = re.compile(name_subfield_regexp)
        subfield_names = [x[0] for x in name_subfield]
        subfield_positions = [int(x[1]) for x in name_subfield]
    else:
        r = None
        subfield_names = None
        subfield_positions = None
    return r, subfield_names, subfield_positions

def name_parse(file, subfield_regexp, subfield_positions):
    if subfield_positions:
        values = []
        s_ = subfield_regexp.search(file)
        for position in subfield_positions:
            values.append(s_.group(position))
        values = tuple(values)
    else:
        values = None
    return values


def name_parse_argumentize(subfield_names, subfield_positions, values):
    return {subfield_names[i-1]:values[i-1] for i in subfield_positions}

########################################################################################################################

class PercentReporter(object):
    def __init__(self, level, total, increment=10, pattern="%i %% complete"):
        self.level = level
        self.total = total
        self.increment = increment
        self.last_reported = 0
        self.pattern = pattern

    def update(self, i, text=None, force=False):
        percent = int(i*100.0/self.total)

        if force or  percent >= self.last_reported + self.increment:
            self.last_reported = percent

            if not text:
                text = self.pattern

            logging.log(self.level, text, percent)

ERROR_REGEXP = re.compile('[^0-9a-zA-Z]+')
########################################################################################################################

def load_list(path):
    k = pandas.read_table(path, header=None)
    return k.iloc[:,0]

def to_dataframe(data, columns, fill_na=None, to_numeric=None):
    if len(data):
        data = pandas.DataFrame(data, columns=columns)
        data = data[columns]
    else:
        data =pandas.DataFrame(columns=columns)
    if to_numeric:
        if type(to_numeric) == list:
            for k in to_numeric:
                data[k] = pandas.to_numeric(data[k])
        elif type(to_numeric) == str:
            for k in data.columns_values:
                data[k] = pandas.to_numeric(data[k], errors=to_numeric)
    if fill_na:
        data = data.fillna(fill_na)
    return data

_gz_re=re.compile(".gz$")
def _compression(path):
    return "gzip" if _gz_re.search(path) else None

def save_dataframe(d, path, mode="w", sep="\t", header=True, fill_na=True, quoting=csv.QUOTE_MINIMAL):
    na_rep = "NA" if fill_na else ''
    compression = _compression(path)
    ensure_requisite_folders(path)
    d.to_csv(path, header=header, mode=mode, compression=compression, sep=sep, index=False, na_rep=na_rep, quoting=quoting)

########################################################################################################################

def open_any(path, mode=None):
    if _gz_re.search(path):
        if mode is None:
            mode = "rb"
        elif mode is "w":
            mode = "wb"
        return gzip.open(path, mode)
    else:
        mode = "r" if mode is None else mode
        return open(path, mode)

def _iterate_file(file, skip_first=False):
    if type(file.mode) == str:
        is_binary = "b" in file.mode
    elif "myfileobj" in dir(file):
        is_binary = "b" in file.myfileobj.mode
    else:
        raise ReportableException("Can't infer mode")

    if skip_first:
        file.readline()

    for i, line in enumerate(file):
        if is_binary:
            line = line.decode()
        yield i, line

def iterate_file(path, skip_first=False):
    with open_any(path) as file:
        for i,line in _iterate_file(file, skip_first):
            yield i, line

def write_to_file(path, lines):
    encode = _gz_re.search(path) is not None
    with open_any(path, "w") as out:
        for line in lines:
            if encode:
                line = line.encode()
            out.write(line)

def to_line(comps):
    return "{}\n".format("\t".join(comps))

def lineify(comps_generator):
    for i,comps in comps_generator:
        yield "\t".join(comps)

def write_iterable_to_file(iterable, path, header=None):
    def to_lines(iterable):
        if header:
            yield header
        for v in iterable:
            yield "{}\n".format(v)
    write_to_file(path, to_lines(iterable))

def get_header(path):
    with open_any(path) as f:
        header = f.readline()
    return header
