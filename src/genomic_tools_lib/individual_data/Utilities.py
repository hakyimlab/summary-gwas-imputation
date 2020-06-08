__author__ = "alvaro barbeira"

import logging

import numpy
import pandas
from pyarrow import parquet as pq

from . import Genotype
from . import Study
from ..Exceptions import ReportableException
from ..miscellaneous import Genomics, PandasHelpers

def get_gene_to_row(gene_annotation):
    return {row.gene_name: i for i, row in enumerate(gene_annotation.itertuples())}

def get_gene_annotation(gene_annotation, gene_to_row, gene_id):
    return gene_annotation.iloc[gene_to_row[gene_id]]

def _load_gene_annotation(path):
    logging.info("Loading gene annotation")
    g = pandas.read_table(path)
    logging.info("Processing")
    g = g.rename(columns={"chr":"chromosome", "start_location":"start", "end_location": "end"})
    if "chr" in g.chromosome.values[0]:
        logging.info("Discarding non-autosomes")
        g = g[g.chromosome.str.contains("chr(\d+)")]
        g = g.assign(chromosome = g.chromosome.str.extract('chr(\d+)').astype(int))
    logging.info("Discarding sex chromosome data, if any")
    g['chromosome'] = pandas.to_numeric(g['chromosome'], errors='coerce')
    g = g.dropna(subset=['chromosome'])
    g.chromosome = g.chromosome.astype(int)
    return g

def _filter_gene_annotation(g, chromosome=None, sub_batches=None, sub_batch=None):
    if chromosome:
        g = g.loc[g.chromosome == chromosome]
    if sub_batches is not None and sub_batch is not None:
        g = PandasHelpers.sub_batch(g, sub_batches, sub_batch)
    return g

def load_gene_annotation(path, chromosome=None, sub_batches=None, sub_batch=None):
    g = _load_gene_annotation(path)
    return _filter_gene_annotation(g, chromosome, sub_batches, sub_batch)

def trim_variant_metadata_on_gene_annotation(vm_, gene_annotation, window, _log_level_v=logging.INFO, id_column="id"):
    variants=set()
    for i, row in enumerate(gene_annotation.itertuples()):
        m_ = Genomics.entries_for_gene_annotation(row, window, vm_)
        if len(m_) == 0:
            logging.log(_log_level_v, "No variants for gene: %s", row.gene_id)
            continue
        variants.update(m_[id_column])

    total = vm_.shape[0]
    vm_ = vm_.loc[vm_[id_column].isin(variants)]
    vm_ = vm_.reset_index(drop=True)
    logging.info("Retained metadata for %d/%d variants", vm_.shape[0], total)
    return vm_

def trim_variants_metadata_on_chromosome(vm_, chromosome):
    vm_ = vm_.loc[vm_.chromosome == chromosome]
    return vm_

def trim_gene_annotation(gene_annotation, variants, window):
    cleared=set()
    for i, row in enumerate(gene_annotation.itertuples()):
        v_ = Genomics.entries_for_gene_annotation(row, window, variants)
        if v_.shape[0] > 0:
            cleared.add(row.gene_id)
    gene_annotation = gene_annotation[gene_annotation.gene_is.isin(cleared)].reset_index(drop=True)
    return gene_annotation

def available_genes(study, gene_annotation):
    g_ = {x for x in study.get_available_pheno_list()}
    return [x for x in gene_annotation.gene_id if x in g_]

def _get_variants_for_gene(gene_id, gene_annotation, gene_to_row, window, study):
    snp_metadata = study.get_variants_metadata()
    g_annotation = get_gene_annotation(gene_annotation, gene_to_row, gene_id)
    variants_metadata = Genomics.entries_for_gene_annotation(g_annotation, window, snp_metadata)
    if not variants_metadata.shape[0]:
        raise ReportableException("no_variants_available")

    variants = study.get_variants(variants_metadata.id.tolist(), to_pandas=True)
    return variants, variants_metadata, g_annotation

def _get_sub_study_for_gene(gene_id, gene_name, gene_annotation, gene_to_row, window, study, rename_pheno="pheno"):
    variants, variants_metadata, gene_annotation = _get_variants_for_gene(gene_name, gene_annotation, gene_to_row, window, study)
    logging.log(5, "Gene Name: {}, variants: {}".format(gene_name, variants.shape[1]))

    phenotype = study.get_phenos([gene_id], to_pandas=True)
    if not gene_id in phenotype:
        raise ReportableException("pheno_not_available")

    covariates=study.get_covariates(to_pandas=True)

    key = gene_id
    if rename_pheno:
        phenotype = phenotype.rename(columns={gene_id:rename_pheno})
        key = rename_pheno

    individuals = study.get_individuals()
    m_ = phenotype.merge(variants, on="individual")
    phenotype = m_[[key]]
    variants = [m_[x] for x in variants_metadata.id]
    genotype = Genotype.Genotype(variants, variants_metadata)
    _study = Study.Study(genotype, phenotype, individuals, covariates)

    return _study

def _get_study_for_gene(context, gene_id, gene_name, rename_pheno="pheno"):
    study = _get_sub_study_for_gene(gene_id,
                                    gene_name,
                                    context.get_gene_annotation(),
                                    context.get_gene_to_row(),
                                    context.get_window(),
                                    context.get_study(),
                                    rename_pheno=rename_pheno)
    return study


def _maf_filter_min_threshold(dosage, metadata, individual_ids, maf_threshold):
    f = metadata[Genotype.MetadataTF.ALLELE_1_FREQUENCY]
    if f>0.5: f = 1-f
    return f<maf_threshold

def _biallelic_filter(dosage, metadata, individual_ids):
    return len(metadata[Genotype.MetadataTF.ALLELE_0]) != 1 or len(metadata[Genotype.MetadataTF.ALLELE_1]) != 1

def is_biallelic_variant(variant):
    comps = variant.split("_")
    if len(comps[2]) != 1: return True
    if len(comps[3]) != 1: return True
    return False

def impute_to_mean_conversion(dosage_raw):
    mean = float(numpy.mean(numpy.array([x for x in dosage_raw if x != "NA"], dtype=numpy.float32)))
    comps = [float(x) if x != "NA" else mean for x in dosage_raw]
    return comps

def gene_annotation_from_parquet_phenotype(pq_fp):
    """

    :param pq_fp: Parquet phenotype file.
    :param region_fp: Text file with LD-independent regions. Columns 'chr', 'start', 'stop'
    :return:
    """
    pheno_columns = pq.ParquetFile(pq_fp).metadata.schema.names
    pheno_columns.remove("individual")
    ll = len(pheno_columns)
    logging.info("Creating gene annotation from {} phenotypes".format(ll))
    df = pandas.DataFrame({'gene_name': pheno_columns,
                               'gene_id': pheno_columns,
                               'gene_type': ['NA'] * ll })

    return df

def expand_gene_annotation_from_regions(df, region_fp, chromosome):
    pheno_columns = df.gene_name.values

    region_df = pandas.read_csv(region_fp, engine='python', sep="\s+")
    region_df['chr'] = region_df['chr'].str.lstrip('chr').astype(int)
    region_df = region_df.rename(mapper={'chr': 'chromosome',
                                         'stop': 'end'},
                                 axis=1)
    region_df['region_id'] = region_df.index
    if chromosome is not None:
        region_df = region_df[region_df['chromosome'] == chromosome]
    logging.info("Expanding gene annotation with {} LD-independent regions".format(region_df.shape[0]))
    pheno_region_lst = []
    for r in region_df.itertuples():
        r_dd = {'region_id': r.region_id,
                'start': r.start,
                'end': r.end,
                'chromosome': r.chromosome}
        region_ = str(r.region_id)
        for pheno in pheno_columns:
            dd = r_dd.copy()
            dd['gene_name'] = pheno + "_" + region_
            dd['gene_id'] = pheno
            dd['gene_type'] = "NA"
            pheno_region_lst.append(dd)
    return pandas.DataFrame(pheno_region_lst)

class _StudyBasedContext:
    def get_study(self): raise RuntimeError("Not implemented")
    def get_available_genes(self): raise  RuntimeError("Not implemented")
    def get_gene_to_row(self): raise RuntimeError("Not implemented")
    def get_gene_annotation(self, gene_id): raise RuntimeError("Not implemented")

    def get_variants_metadata(self): raise RuntimeError("Not implemented")
    def get_variants(self, variants): raise RuntimeError("Not implemented")
    def get_variants_metadata_for_gene(self, gene_id): raise RuntimeError("Not implemented")
    def get_variants_for_gene(self, gene_id): raise RuntimeError("Not implemented")
    def get_window(self): raise RuntimeError("Not implemented")


class StudyBasedContext(_StudyBasedContext):
    def __init__(self, study, gene_annotation, window):
        self.study = study
        self.gene_annotation = gene_annotation
        self.gene_to_row = get_gene_to_row(gene_annotation) if gene_annotation is not None else None
        self.window = window
        self.available_genes_ = None

    def get_study(self): return self.study
    def get_window(self): return self.window
    def get_gene_to_row(self): return self.gene_to_row

    def get_gene_annotation(self, gene_id=None):
        if gene_id: return get_gene_annotation(self.gene_annotation, self.gene_to_row, gene_id)
        return self.gene_annotation

    def get_available_genes(self):
        if not self.available_genes_: self.available_genes_ = available_genes(self.study, self.gene_annotation)
        return self.available_genes_

    def get_variants_metadata_for_gene(self, gene_id):
        a = self.get_gene_annotation(gene_id)
        return Genomics.entries_for_gene_annotation(a, self.window, self.study.get_variants_metadata())

    def get_variants_metadata(self):
        return self.study.get_variants_metadata()

    def get_variants(self, variants):
        return self.study.get_variants(variants)