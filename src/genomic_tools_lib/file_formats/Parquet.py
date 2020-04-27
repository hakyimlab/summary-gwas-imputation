__author__ = "alvaro barbeira"

import logging
import pandas
import pyarrow as pa
import pyarrow.parquet as pq
import numpy
import re

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
        return [individuals.index(x) for x in specific_individuals]
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
            indexes = { str(i) for i in specific_individuals}
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

###############################################################################
class MultiFileGenoHandler:
    """
    This class is for loading parquet metadata and genotype files. Most of its
    functionality is meant to assist in the case that both metadata and genotype
    are split into 22 files, but it should be robust to the case where there is
    only one genotype or metadata file.
    """
    def __init__(self, features, metadata):
        """
        If either argument is a pattern for multiple files, it must be
        formattable with the argument 'chr'

        :param features: filepath (or pattern for multiple) for genotype files
        :param metadata: filepath (or pattern for multiple) for geno metadata
        """
        self.m_features = self.check_if_formattable(features)
        if self.m_features:
            self.features = self.format_chrom_file_names(features)
        else:
            self.features = [features]

        self.m_metadata = self.check_if_formattable(metadata)
        if self.m_metadata:
            self.metadata = self.format_chrom_file_names(metadata)
        else:
            self.metadata = [metadata]

    @staticmethod
    def format_chrom_file_names(s):
        l = [s.format(chr=i) for i in range(1, 23)]
        return l
    @staticmethod
    def check_if_formattable(s):
        matches = re.findall('{(.*?)}', s)
        if len(matches) > 0 and matches[0] == 'chr':
            return True
        else:
            return False

    def load_metadata(self, whitelist=None):
        df_lst = []
        for i in self.metadata:
            df_i = pq.read_table(i).to_pandas()
            if whitelist is not None:
                df_i = df_i[df_i.id.isin(whitelist)]
            df_lst.append(df_i)
        return pandas.concat(df_lst)

    def load_features(self, metadata, individuals, pandas=False):
        """
        :param metadata: pandas DataFrame with columns 'variant_id' and
                'chromosome'
        :param individuals: list. Individual IDs
        :param pandas: bool. default False. Whether the returned obj should be
                a pandas DataFrame
        :return: dict.
        """
        if self.m_features:
            d_ = self._load_features_multiple(metadata, individuals, pandas)
        else:
            d_ = self._load_features_single(metadata, individuals, pandas)
        n_requested = metadata.shape[0]
        n_loaded = len(d_) - 1
        if n_requested == n_loaded:
            logging.log(5, "Loaded {} features".format(n_loaded))
            return (d_, metadata)
        else:
            s = "Requested {} features and loaded {}".format(n_requested,
                                                             n_loaded)
            logging.warning(s)
            if pandas:
                loaded_features = set(d_.columns)
            else:
                loaded_features = {v for v in d_.keys()}
            loaded_features.remove('individual')
            loaded_metadata = metadata[metadata.id.isin(loaded_features)]
            return (d_, loaded_metadata)

    def _load_features_single(self, metadata, individuals, pandas):
        dd =  _read(pq.ParquetFile(self.features[0]),
                    columns=[x for x in metadata.id],
                    specific_individuals=individuals,
                    to_pandas=pandas)
        logging.log(5, "Loaded {} features".format(len(dd) - 1))
        return dd

    def _load_features_multiple(self, metadata, individuals, pandas):
        df_lst = []
        i_ = list(individuals)
        for chr, group in metadata.groupby('chromosome'):
            chr_fp = self.features[chr - 1]
            chr_vars = list(group.id)
            chr_features = _read(pq.ParquetFile(chr_fp), chr_vars,
                                         specific_individuals=i_,
                                         to_pandas = True)
            df_lst.append(chr_features.set_index('individual'))
        while len(df_lst) > 1:
            df_lst[0].join(df_lst.pop(), how='inner')
        logging.log(5, "Loaded {} features".format(df_lst[0].shape[1]))
        if pandas:
            return df_lst[0].reset_index()
        else:
            return df_lst[0].reset_index().to_dict(orient='list')

###############################################################################

class PhenoDataHandler:
    """
    This class is meant to handle phenotype data and optional weights.
    """
    def __init__(self, data_fp, data_annotation=None, sub_batches=None,
                 sub_batch=None):
        """

        :param data_fp: str. filepath to a Parquet phenotype file
        :param data_annotation: pandas DataFrame. with columns 'gene_name',
                    'gene_id', 'gene_type'
        :param sub_batches: int. How many batches the data should be split into
        :param sub_batch: int. 0-indexed which batch should be considered
        """
        self.data = pq.ParquetFile(data_fp)
        d_names = self.data.metadata.schema.names
        d_names.remove('individual')
        if data_annotation is not None:
            self.data_annotation = data_annotation
        else:
            self.data_annotation = self._load_da_manual(sub_batches,
                                                        sub_batch)

        self._features_metadata = None
        self.send_weights = False
        self._features_weights = None
        self.features = None

    def _load_da_manual(self, n_batches = None, batch = None):
        names = self.data.metadata.schema.names
        if (n_batches is not None) and (batch is not None):
            names = numpy.array_split(names, n_batches)[batch]

        annot_dd = {'gene_name': names, 'gene_id': names,
                    'gene_type': ['NA'] * len(names)}
        df = pandas.DataFrame(annot_dd)
        return df

    def add_features_metadata(self, metadata):
        """

        :param metadata: pandas DataFrame
        """
        self._features_metadata = metadata
        self._merge_metadata_weights()

    def add_features_weights(self, weights, pheno_col='gene_id'):
        """
        :param weights: pandas DataFrame
        :param pheno_col: str. Name of col with phenotype data
        """

        weights = weights.rename(mapper={'variant_id': 'id', pheno_col: 'gene_id'}, axis=1)

        # Only keep the intersection of genes from weights, genes from data
        ids_from_weights = set(weights.gene_id)
        ids_from_data = set(self.data_annotation.gene_id)
        gene_ids = ids_from_weights.intersection(ids_from_data)
        weights = weights.loc[weights.gene_id.isin(gene_ids)]
        self.data_annotation = \
            self.data_annotation.loc[self.data_annotation.gene_id.isin(gene_ids)]
        logging.info("Weights loaded for {} genes".format(len(gene_ids)))

        self._features_weights = weights
        self._merge_metadata_weights()

    def _merge_metadata_weights(self):
        if (self._features_weights is not None) and \
        (self._features_metadata is not None):
            logging.info( 'Merging geno metadata and weights')
            f_w = self._features_weights.set_index('id')
            m = self._features_metadata.set_index('id')
            cols = list(m.columns)
            cols.extend(['w', 'gene_id'])
            f_w = f_w.join(m, how='left', rsuffix='_m')[cols]
            f_w['id'] = f_w.index
            self._features_weights = f_w.groupby('gene_id')
            self.send_weights = True

    def get_features(self, pheno):
        if self.send_weights:
            return self._features_weights.get_group(pheno)
        else:
            return self._features_metadata[['id', 'chromosome']]

    def load_pheno(self, pheno, to_pandas=False):
        return _read(self.data, [pheno], to_pandas=to_pandas)

    def get_individuals(self):
        return self.data.read(columns=['individual']).to_pydict()['individual']
