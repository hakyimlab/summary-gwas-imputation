__author__ = "alvaro barbeira"

import logging
import pandas
import numpy

from . import SummaryInputation
from ..data_management import TextFileTools
from ..file_formats import Parquet
from ..file_formats.gwas import Utilities as GWASUtilities
from ..individual_data.Utilities import StudyBasedContext
from ..miscellaneous import Genomics
from ..miscellaneous import PandasHelpers

########################################################################################################################
class _Mixin:
    def __init__(self, gwas, cutoff, regularization, frequency_filter, standardise_dosages, keep_palindromic_imputation, use_palindromic_snps):
        self.gwas_data = _parse_gwas(gwas)
        self.cutoff = cutoff
        self.regularization = regularization
        self.frequency_filter = frequency_filter
        self._standardise_dosages = standardise_dosages
        self._keep_palindromic_imputation = keep_palindromic_imputation
        self._use_palindromic_snps = use_palindromic_snps

    def get_gwas_slice(self, variants_metadata=None, variant=None):
        return _gwas_for_slice(self.gwas_data, variants_metadata, variant)

    def get_cutoff(self): return self.cutoff
    def get_regularization(self): return self.regularization
    def get_freq_filter(self): return self.frequency_filter
    def standardise_dosages(self): return self._standardise_dosages
    def keep_palindromic_imputation(self): return self._keep_palindromic_imputation
    def use_palindromic_snps(self): return self._use_palindromic_snps

class VariantContext(StudyBasedContext, _Mixin, SummaryInputation._VariantContext):
    def __init__(self, study, window, gwas, cutoff, regularization, frequency_filter, standardise_dosages, specific_target_variants=None, keep_palindromic_imputation=None, use_palindromic_snps=None):
        super().__init__(study, None, window)
        _Mixin.__init__(self, gwas, cutoff, regularization, frequency_filter, standardise_dosages, keep_palindromic_imputation, use_palindromic_snps)
        self.target_variants_metadata = specific_target_variants if specific_target_variants is not None else study.get_variants_metadata()

    def get_target_variants_metadata(self):
        return self.target_variants_metadata

class CachingVariantContext(VariantContext):
    def __init__(self, study, window, gwas, cutoff, regularization, frequency_filter, standardise_dosages, specific_target_variants=None, keep_palindromic_imputation=None, use_palindromic_snps=None):
        super().__init__(study, window, gwas, cutoff, regularization, frequency_filter, standardise_dosages, specific_target_variants, keep_palindromic_imputation, use_palindromic_snps)
        self.cache = {}

    def get_variants(self, variants):
        return _get_variants_cached(variants, self.study, self.cache)

def _get_variants_cached(variants, study, cache):
    req = [x for x in variants if not x in cache]
    if len(req):
        v = study.get_variants(req)
        for k in v:
            cache[k] = v[k].values
    return {k:cache[k] for k in variants}

########################################################################################################################

class RegionContext(StudyBasedContext, _Mixin, SummaryInputation._RegionContext):
    def __init__(self, study, window, gwas, cutoff, regularization, frequency_filter, regions, standardise_dosages=None, keep_palindromic_imputation=None, use_palindromic_snps=None):
        super().__init__(study, None, window)
        _Mixin.__init__(self, gwas, cutoff, regularization, frequency_filter, standardise_dosages, keep_palindromic_imputation, use_palindromic_snps)
        self.target_regions = regions

    def get_target_regions(self):
        return self.target_regions


########################################################################################################################

def _parse_gwas(gwas):
    gwas = gwas[~gwas.zscore.isna()]
    logging.log(9, "Parsing GWAS")
    g = {}
    for i,t in enumerate(gwas.itertuples()):
        if t.zscore == "NA":
            continue

        if not t.chromosome in g:
            g[t.chromosome] = {}

        c = g[t.chromosome]

        #Yeah, here is where pandas' reading the position column as a float bites us
        if numpy.isnan(t.position): continue
        _pos = int(t.position)
        if not _pos in c:
            c[_pos] = {}
        c = c[_pos]

        z = float(t.zscore)
        c[t.effect_allele] = {}
        c[t.effect_allele][t.non_effect_allele] = (z, t.variant_id)

        c[t.non_effect_allele] = {}
        c[t.non_effect_allele][t.effect_allele] = (-z, t.variant_id)
    return g

def _gwas_for_slice(gwas_data, variants_metadata, variant):
    if variants_metadata is None and variant is not None:
        variants_metadata = pandas.DataFrame(
            [(variant.id, variant.chromosome, variant.position, variant.non_effect_allele, variant.effect_allele, variant.effect_allele_frequency)],
            columns=["id", "chromosome", "position", "non_effect_allele", "effect_allele", "effect_allele_frequency"])

    d = []
    for i,t in enumerate(variants_metadata.itertuples()):
        chr = "chr{}".format(t.chromosome)

        if not chr in gwas_data: continue
        c = gwas_data[chr]

        if not t.position in c: continue
        c = c[t.position]

        if not t.effect_allele in c: continue
        c = c[t.effect_allele]

        if not t.non_effect_allele in c: continue
        zscore, id = c[t.non_effect_allele]
        #Use the GTEx id? id = t.rsid

        d.append((id, t.id, t.chromosome, t.position, t.effect_allele, t.non_effect_allele, t.effect_allele_frequency, zscore, None, "measured"))
    return pandas.DataFrame(d, columns=SummaryInputation.Results._fields)

def load_gwas(args, use_specific_targets):
    """Watch out! Pandas parser on occasion reads the column `position` as integer. Yikes."""
    columns = ["variant_id", "panel_variant_id", "chromosome", "position", "non_effect_allele", "effect_allele","zscore"]
    if use_specific_targets is not None:
        logging.info("Acquiring filter tree for %d targets", use_specific_targets.shape[0])
        tree = GWASUtilities.get_chromosome_position_tree(use_specific_targets)
        filter = GWASUtilities.get_filter(tree)
        logging.info("Processing gwas source")

        gwas = TextFileTools.load_dataframe(args.gwas_file, additional_filter=filter, columns=columns)[columns]
    else:
        logging.info("opening gwas file")
        gwas = pandas.read_table(args.gwas_file, usecols =columns)[columns]
        if args.chromosome: gwas = gwas.loc[gwas.chromosome == "chr{}".format(args.chromosome)]

    logging.log(9, "Loaded %d GWAS variants", gwas.shape[0])
    #Mind, the imputation must do the discarding of ambiguous
    #gwas = GWASUtilities.discard_ambiguous(gwas)

    return  gwas

def trim_variant_metadata(args, vm, use_specific_targets):
    chromosomes = {x for x in use_specific_targets.chromosome}
    r = []
    for chromosome in chromosomes:
        u = use_specific_targets[use_specific_targets.chromosome == chromosome]
        window_start = u.position.min() - args.window
        window_end = u.position.max() + args.window
        v = Genomics.entries_for_window(chromosome, window_start, window_end, vm)
        r.append(v)
    return pandas.concat(r)

def load_study(args):
    study = Parquet.study_from_parquet(args.parquet_genotype, args.parquet_genotype_metadata, chromosome=args.chromosome)
    vm = study.variant_metadata
    vm.rename(columns={"allele_0":"non_effect_allele", "allele_1":"effect_allele", "allele_1_frequency":"effect_allele_frequency"}, inplace=True)
    study.variant_metadata = vm[['chromosome', 'position', 'id', 'non_effect_allele', 'effect_allele', 'effect_allele_frequency', 'rsid']]
    return study

########################################################################################################################

def context_from_args(args):
    logging.info("Creating context by variant")

    logging.info("Loading study")
    study = load_study(args)

    if args.sub_batches and args.sub_batch is not None:
        logging.log(9, "Selecting targets from sub-batches")
        use_specific_targets = PandasHelpers.sub_batch(study.variant_metadata, args.sub_batches, args.sub_batch)
        study.variant_metadata = trim_variant_metadata(args, study.variant_metadata, use_specific_targets)
    else:
        use_specific_targets = None

    logging.info("Loading gwas")
    gwas = load_gwas(args,  use_specific_targets)

    if args.cache_variants:
        context = CachingVariantContext(study, args.window, gwas, args.cutoff, args.regularization, args.frequency_filter, args.standardise_dosages, use_specific_targets)
    else:
        context = VariantContext(study, args.window, gwas, args.cutoff, args.regularization, args.frequency_filter, args.standardise_dosages, use_specific_targets)
    return context

########################################################################################################################

def load_region(args, study):
    regions = pandas.read_table(args.by_region_file).rename(columns={"chr":"chromosome", "stop":"end"})
    regions = regions.dropna()
    regions = regions.assign(chromosome = regions.chromosome.str.split("chr").str.get(1).astype(numpy.int32))

    use_specific_targets = None

    if args.chromosome:
        logging.log(9, "Selecting target regions with specific chromosome")
        regions = regions.loc[regions.chromosome == args.chromosome]

        if args.containing:
            logging.log(9, "Selecting target regions with specific position")
            regions = regions[(regions.chromosome == args.chromosome) & (regions.start <= args.containing) & (args.containing <= regions.end)]
            regions = regions.reset_index(drop=True)

    if args.sub_batches and args.sub_batch is not None:
        logging.log(9, "Selecting target regions from sub-batches")
        regions = PandasHelpers.sub_batch(regions, args.sub_batches, args.sub_batch)

    if args.chromosome or (args.sub_batches and args.sub_batch is not None):
        logging.log(9, "generating GWAS whitelist")
        use_specific_targets = []
        for region in regions.itertuples():
            w = Genomics.entries_for_window(region.chromosome, region.start - args.window, region.end + args.window, study.variant_metadata)
            use_specific_targets.append(w)
        use_specific_targets = pandas.concat(use_specific_targets).drop_duplicates()

    return regions, use_specific_targets


def context_by_region_from_args(args):
    logging.info("Creating context by variant")

    logging.info("Loading study")
    study = load_study(args)

    logging.info("Loading regions")
    regions, use_specific_targets = load_region(args, study)

    logging.info("Loading gwas")
    gwas = load_gwas(args, use_specific_targets)
    context = RegionContext(study, args.window, gwas, args.cutoff, args.regularization, args.frequency_filter, regions,
                            args.standardise_dosages, args.keep_palindromic_imputation, args.use_palindromic_snps)

    return context

