import gzip
import pandas
import logging

from .TextFileTools import sanitize_components
from .. import Utilities

def data_frame_streamer(path, sentinel_column, sentinel_white_list=None, sentinel_suffix=None, additional_skip_row_check=None, to_numeric=None, sanitize=False, header=None):
    _check_sentinel(sentinel_white_list, sentinel_suffix)

    found = set()

    sentinel_index = None
    sentinel = None
    buffer = []
    for i, line in Utilities.iterate_file(path):
        comps = line.strip().split()

        if i == 0:
            if not header:
                header = comps
                sentinel_index = header.index(sentinel_column)
                continue
            else:
                sentinel_index = header.index(sentinel_column)

        comps = tuple(comps)

        row_sentinel = comps[sentinel_index]

        if _should_skip(row_sentinel, sentinel_white_list, sentinel_suffix):
            continue

        if additional_skip_row_check and additional_skip_row_check(comps):
            continue

        if not sentinel:
            sentinel = row_sentinel

        if sanitize:
            comps = sanitize_components(comps)

        if row_sentinel == sentinel:
            buffer.append(comps)
            continue

        data = Utilities.to_dataframe(buffer, header, to_numeric=to_numeric)
        yield data

        if sentinel_white_list is not None:
            found.add(sentinel)
            if len(found) == len(sentinel_white_list):
                logging.info("Found everything in the whitelist")
                return

        buffer = [comps]
        sentinel = comps[sentinel_index]

    if len(buffer):
        data = Utilities.to_dataframe(buffer, header, to_numeric=to_numeric)
        yield data


def _check_sentinel(sentinel_white_list, sentinel_suffix):
    if not sentinel_white_list:
        return

    if sentinel_suffix:
        w = [x for x in sentinel_white_list if sentinel_suffix in x]
        if len(w):
            raise RuntimeError("Unexpected condition: sentinel_whitelist values contain suffix. Verify the white list, or don't use the suffix")

def _should_skip(row_sentinel, sentinel_white_list, sentinel_suffix):
    if sentinel_white_list is not None:
        ward = row_sentinel
        if sentinel_suffix:
            ward = ward.split(sentinel_suffix)[0]

        if not ward in sentinel_white_list:
            return True

    return False

def load_filtered_data_frame(path, sentinel_column, sentinel_white_list=None, sentinel_suffix=None, rename_columns=None, columns_filter=None):
    result = pandas.DataFrame()
    for i, d in enumerate(data_frame_streamer(path, sentinel_column, sentinel_white_list, sentinel_suffix)):
        if rename_columns:
            d = d.rename(columns=rename_columns)
        if columns_filter:
            d = d[columns_filter]
        result = pandas.concat([result, d])
    return result