__author__ = "alvaro barbeira"

import logging
import pandas
import pyarrow as pa
import pyarrow.parquet as pq
import numpy

from ..individual_data import Study, Genotype
from ..DataSink import DataFrameSink
from ..miscellaneous import Genomics

def _deplete_variants_to_record_batch(variants, variant_ids, individual_ids):
    data = [pa.array(individual_ids)]
    logging.log(8, "Consuming genotype variant data")
    for i in range(0, len(variants)):
        variant = variants.pop(0)
        data.append(pa.array(variant))
    logging.log(8, "Consumed genotype variant data")
    names =  ["individual"]+variant_ids
    return pa.RecordBatch.from_arrays(data, names)

def _deplete_genotype_variants_to_record_batch(genotype, individual_ids):
    return _deplete_variants_to_record_batch(genotype.variants, [x for x in genotype.metadata.id.values], individual_ids)

#The conversion to record batch is because of a bug in pyarrow for flavor="spark"
def _to_record_batch(df):
    data =[]
    names = list(df.columns.values)
    for c in names:
        data.append(pa.array(df[c]))
    return  pa.RecordBatch.from_arrays(data, names)

def save_variants(path, genotype, individual_ids):
    """Will consume the data in the study"""
    record_batch = _deplete_genotype_variants_to_record_batch(genotype, individual_ids)
    table = pa.Table.from_batches([record_batch])
    pq.write_table(table, path, flavor="spark")

def _save_metadata(path, metadata):
    table = _to_record_batch(metadata.iloc[0:2, ])
    with ParquetDataFrameSink(path, table.schema) as sink:
        for c_ in range(1, 23):
            logging.log(8, "Saving metadata for chromosome %d", c_)
            p_ = metadata.loc[metadata.chromosome == c_]
            sink.sink(p_)

def save_metadata(path, genotype):
    _save_metadata(path, genotype.metadata)

def save_variable(path, variables, individual_ids=None):
    if individual_ids:
        variables = variables.copy()
        columns = list(variables.columns.values)
        columns = ["individual"] + columns
        variables["individual"] = individual_ids
        variables = variables[columns]
    batches = [_to_record_batch(variables)]
    table = pa.Table.from_batches(batches)
    pq.write_table(table, path, flavor="spark")

def save_study(study, path_prefix):
    """Will consume the data in the study"""
    genotype = study.get_genotype()
    individuals = study.get_individuals()
    phenotype = study.get_phenotype()
    covariates = study.get_covariates()

    path_variant = path_prefix + ".variants.parquet"
    save_variants(path_variant, genotype, individuals)

    path_variant_metadata = path_prefix + ".variants_metadata.parquet"
    save_metadata(path_variant_metadata, genotype)

    path_pheno = path_prefix + ".pheno.parquet"
    save_variable(path_pheno, phenotype, individuals)

    if covariates is not None:
        path_covariate = path_prefix + ".covariate.parquet"
        save_variable(path_covariate, covariates, individuals)

########################################################################################################################
def variant_key_value_from_metadata(path):
    p = pq.read_table(path, columns=["id", "rsid"]).to_pandas()
    return {x.id:x.rsid for x in p.itertuples()}

def variants_from_metadata(path, frequency_threshold=None):
    if not frequency_threshold:
        p = pq.read_table(path, columns=["id"]).to_pandas()
    else:
        p = pq.read_table(path, columns=["id", "allele_1_frequency"]).to_pandas().rename(columns={"allele_1_frequency":"f"})
        p = p.loc[(p.f > frequency_threshold) & (p.f < (1-frequency_threshold))]
    return {x for x in p.id.values}

########################################################################################################################

class ParquetStudyBase(Study._Study):
    def __init__(self, variant_metadata, phenotype_file, covariates, individuals):
        super().__init__()
        self.variant_metadata = variant_metadata
        self.phenotype_file = phenotype_file
        self.individuals = individuals
        self.pheno_list_ = None
        self.covariates = covariates

    def get_variants_metadata(self, variants=None): return Genotype._get_variants_metadata(self.variant_metadata, variants)
    def get_phenos(self, phenos=None, to_pandas=True): return _read(self.phenotype_file, phenos, to_pandas=to_pandas)
    def get_individuals(self): return self.individuals

    #def get_available_pheno_list(self): return Study._get_list(self.pheno)
    def get_available_covariate_list(self): return Study._get_list(self.covariates)
    def get_covariates(self, covariates=None, to_pandas=True): return Study._get(self.covariates, covariates, to_pandas=to_pandas)

    def get_available_pheno_list(self):
        if not self.pheno_list_:
            self.pheno_list_ = [x for x in self.phenotype_file.schema.names if x != "individual"]
        return self.pheno_list_

class ParquetStudy(ParquetStudyBase):
    """
variant_file and phenotype_file are meant to be read on the fly.
Variant data will be loaded into a dictionary to avoid pandas overhead. This is the opposite default to other sibling Study classes.
variant_metadata, we'll preload.
    """
    def __init__(self, variant_file, variant_metadata, phenotype_file, covariates, individuals):
        super().__init__(variant_metadata, phenotype_file, covariates, individuals)
        self.variant_file = variant_file

    def get_variants(self, variants=None, to_pandas=False, omit_individuals=False, specific_individuals=None): return _read(self.variant_file, variants, omit_individuals, to_pandas, specific_individuals)

########################################################################################################################

class ParquetSplitStudy(ParquetStudyBase):
    """
variant_file_map is a map from chromosome numbers to parquet variant filed
phenotype_file is meant to be read on the fly.
Variant data will be loaded into a dictionary to avoid pandas overhead. This is the opposite default to other sibling Study classes.
variant_metadata, we'll preload.
    """
    def __init__(self, variant_file_map, variant_metadata, phenotype_file=None, covariates=None, individuals=None):
        super().__init__(variant_metadata, phenotype_file, covariates, individuals)
        self.variant_file_map = variant_file_map

    def get_variants(self, variants=None, to_pandas=False, omit_individuals=False, specific_individuals=None):
        """Asssumes all requested varianst in a same chromosome"""
        if variants is None:
            raise RuntimeError("This implementation demands a specific list of variants")
        chr = variants[0].split("_")[0]
        v = self.variant_file_map[chr]
        return _read(v, variants, omit_individuals, to_pandas, specific_individuals)

class ParquetSingleSplitStudy(ParquetStudyBase):
    """
variant_file_map is a map from chromosome numbers to paths. Each genotype file will be opened when a variant requirement needs it.
Support only one chromosome open at any givven time.
phenotype_file is meant to be read on the fly.
Variant data will be loaded into a dictionary to avoid pandas overhead. This is the opposite default to other sibling Study classes.
variant_metadata, we'll preload.
    """
    def __init__(self, variant_paths, variant_metadata, phenotype_file=None, covariates=None, individuals=None):
        super().__init__(variant_metadata, phenotype_file, covariates, individuals)
        self.variant_paths = variant_paths
        self.last_chr = None
        self.file = None

    def get_variants(self, variants=None, to_pandas=False, omit_individuals=False, specific_individuals=None):
        """Asssumes all requested varianst in a same chromosome"""
        if variants is None:
            raise RuntimeError("This implementation demands a specific list of variants")
        chr = variants[0].split("_")[0]
        if chr != self.last_chr:
            logging.log(9, "Loading new chromosome requirement: %s", chr)
            self.last_chr = chr
            path = self.variant_paths[chr]
            self.file = pq.ParquetFile(path)
            logging.log(9, "Loaded new chromosome requirement: %s", chr)
        return _read(self.file, variants, omit_individuals, to_pandas, specific_individuals)

########################################################################################################################

def _individual_mask(individuals, specific_individuals):
    if specific_individuals:
        available_specific_individuals = len([x for x in specific_individuals if x in individuals])
        if available_specific_individuals < len(specific_individuals):
            logging.warning("Requested specific individuals are missing from the genotype")
        return [individuals.index(x) for x in specific_individuals if x in individuals]
    else:
        return None

def get_snps_data(annotation_row, window, snp_metadata, snp_file, specific_individuals, to_pandas=False):
    features_in_window = Genomics.entries_for_gene_annotation(annotation_row, window, snp_metadata)
    return features_in_window, _read(snp_file, [x for x in features_in_window.id.values], to_pandas=to_pandas, specific_individuals=specific_individuals)

def _read(file, columns=None, skip_individuals=False, to_pandas=False, specific_individuals=None):
    if columns is None:
        columns = file.schema.names
    if not skip_individuals:
        columns = ["individual"]+[x for x in columns]
    if skip_individuals and specific_individuals is not None:
            raise RuntimeError("Unsupported combination")
    v = file.read(columns=columns)
    if to_pandas:
        v = v.to_pandas()
        if specific_individuals is not None:
            indexes = set(specific_individuals)
            v = v.loc[v.individual.isin(indexes)]
    else:
        if specific_individuals:
            mask = _individual_mask(v.column(0).to_pylist(), specific_individuals)

        v = {c.name:(numpy.array(c.to_pylist(), dtype=numpy.float32) if c.name != "individual" else numpy.array(c.to_pylist(), dtype=numpy.str)) for c in v}
        if specific_individuals:
            v = {k:d[mask] for k,d in v.items()}
    return v

def study_from_parquet(variants, variants_metadata, pheno=None, covariates=None, post_process_variants_metadata=None, chromosome=None, frequency_filter=None):
    logging.info("Loading variants' parquet file")
    _v = pq.ParquetFile(variants)

    logging.info("Loading variants metadata")
    if post_process_variants_metadata or chromosome:
        f = pq.ParquetFile(variants_metadata)
        _vm =[]
        if chromosome:
            _r = (chromosome-1,)
        else:
            _r = f.num_row_groups if f.num_row_groups < 22 else 22
            _r = range(0, _r)

        for i in _r:
            logging.log(9,"Loading row group %d", i)
            vmi_ = f.read_row_group(i).to_pandas()

            if post_process_variants_metadata:
                vmi_ = post_process_variants_metadata(vmi_)

            if vmi_.shape[0] > 0:
                _vm.append(vmi_)

        _vm = pandas.concat(_vm)
    else:
        _vm = pq.ParquetFile(variants_metadata).read().to_pandas()

    if frequency_filter:
        _vm = _vm.loc[(_vm.allele_1_frequency > frequency_filter) & (_vm.allele_1_frequency < (1 - frequency_filter))]

    if pheno:
        logging.info("Loading phenotype")
        p_ = pq.ParquetFile(pheno)
    else:
        p_ = None
    individuals = _v.read(["individual"]).to_pandas()["individual"]

    if covariates:
        logging.info("Loading covariates")
        c_ = pq.ParquetFile(covariates).read().to_pandas()
    else:
        c_ = None

    return ParquetStudy(_v, _vm, p_, c_, individuals)

########################################################################################################################

class ParquetDataFrameSink(DataFrameSink):
    def __init__(self, path, schema, compression=None):
        self.path = path
        self.schema = schema
        self.writer = None
        self.compression = compression

    def __enter__(self):
        self.initialize()
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.finalize()

    def sink(self, d):
        d = _to_record_batch(d)
        d = pa.Table.from_batches([d])
        self.writer.write_table(d)

    def initialize(self):
        logging.log(9, "Initializing parquet sink")
        self.writer = pq.ParquetWriter(self.path, self.schema, flavor="spark", compression=self.compression)

    def finalize(self):
        logging.log(9, "Finalizing parquet sink")
        self.writer.close()