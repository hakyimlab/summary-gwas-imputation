__author__ = "alvaro barbeira"

from . import Exceptions

class DataSink:
    def sink(self, d):
        raise Exceptions.ReportableException("Not implemented")

    def initialize(self):
        raise Exceptions.ReportableException("Not implemented")

    def finalize(self):
        raise Exceptions.ReportableException("Not implemented")

    def __enter__(self):
        raise Exceptions.ReportableException("Not implemented")

    def __exit__(self, exc_type, exc_val, exc_tb):
        raise Exceptions.ReportableException("Not implemented")

class DataFrameSink(DataSink):
    pass
