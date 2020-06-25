import os
import logging
from timeit import default_timer as timer
import re
import sys

from genomic_tools_lib import Logging, Utilities
from genomic_tools_lib.file_formats import Parquet

from genomic_tools_lib.external_tools.tensorqtl import genotypeio, trans
from pyarrow import parquet as pq
import pandas
import numpy

class FileOut:
    def __init__(self, out_dir, chromosome, parquet_geno_metadata=None):
        self.out_dir = out_dir
        self.chromosome = chromosome
        self.dir_pattern = "chr-{chr}".format(chr=chromosome)
        Utilities.maybe_create_folder(os.path.join(out_dir, self.dir_pattern))
        self.file_pattern = "tensorqtl-summ-stats_{pheno}_chr{chr}.txt"
        self.GWAS_COLS = ['variant_id', 'panel_variant_id', 'chromosome',
                      'position', 'effect_allele', 'non_effect_allele',
                      'current_build', 'frequency', 'sample_size', 'zscore',
                      'pvalue', 'effect_size', 'standard_error',
                      'imputation_status', 'n_cases', 'gene', 'region_id']
        if parquet_geno_metadata is not None:
            self.metadata = pq.read_table(parquet_geno_metadata)
            self.metadata.set_index('variant_id')
        else:
            self.metadata = None
    def write_text(self, df, pheno):
        fp = os.path.join(self.out_dir,
                          self.dir_pattern,
                          self.file_pattern.format(pheno=pheno,
                                                   chr=self.chromosome))
        df.to_csv(fp, sep = "\t", index=False)

    def make_gwas(self, df, n_samples):
        """
        Take tensorqtl output, with cols
        ['variant_id', 'phenotype_id', 'pval', 'b', 'b_se', 'maf']
        and create a GWAS file.
        """
        ll = df.shape[0]
        df.index = df['variant_id']
        map_dd = {'phenotype_id': 'gene',
                  'pval': 'pvalue',
                  'b': 'effect_size',
                  'b_se': 'standard_error',
                  'maf': 'frequency'}
        df = df.rename(mapper=map_dd, axis=1)
        df['sample_size'] = [n_samples] * ll
        if self.metadata is not None:
            df = df.merge(self.metadata, how='left')
        return self._fill_empty_gwas_cols(df)

    def _fill_empty_gwas_cols(self, df, fill='NA'):
        ll = df.shape[0]
        for i in self.GWAS_COLS:
            if i not in df:
                logging.log(7, "Filling column {}".format(i))
                df[i] = [fill] * ll
        return df[self.GWAS_COLS]


class FileIn:
    def __init__(self, plink_geno, parquet_pheno, rewrite_bim=False):
        self.BIM = ".bim"
        self.BED = ".bed"
        self.FAM = ".fam"
        self.geno_pre = plink_geno
        if rewrite_bim:
            self._rewrite_bim()
        self.pheno_fp = parquet_pheno
        self.pheno = pq.ParquetFile(parquet_pheno)

        self.covar_df = None

        self.individuals = self._find_intersection()



    def get_pheno(self): return self._load_pheno(self.individuals)
    def get_covar(self): return self.covar_df

    def get_geno(self):
        logging.info("Loading genotype")
        pr = genotypeio.PlinkReader(self.geno_pre, verbose=False)
        genotype_df = pr.load_genotypes()
        genotype_df.columns = [str(i) for i in genotype_df.columns]
        logging.info("Finished loading genotype")

        return genotype_df

    @staticmethod
    def _load_fam(fp):
        df = pandas.read_csv(fp,
                        sep='\s',
                        header=None,
                        names=['FID', 'IID', 'PID', 'MID', 'SEX', 'PHENO'],
                        engine='python')
        return df

    def _rewrite_bim(self):
        logging.info("Rewriting BIM to have standard variant IDs")
        fp = self.geno_pre + self.BIM
        df = pandas.read_csv(fp,
                             sep="\s",
                             header=None,
                             names=['CHR', 'SNP', 'XXX', 'BP', 'REF', 'ALT'])
        df['SNP'] = ("chr" + df.CHR.astype('str')
                     + "_" + df.BP.astype('str')
                     + "_" + df.REF
                     + "_" + df.ALT)

        df.to_csv(fp, sep="\t")

    def _find_intersection(self):
        fam_df = self._load_fam(self.geno_pre + self.FAM)

        geno_ids = set(fam_df.IID.astype(str).values)
        logging.log(9, "Individuals from genotype: {}".format(len(geno_ids)))

        pheno_ids = set(Parquet._read(self.pheno, columns=['individual'])['individual'])
        logging.log(9, "Individuals from phenotype: {}".format(len(pheno_ids)))

        # I have so much trouble remembering python set methods.
        # This is intersection.
        all_ids = geno_ids.intersection(pheno_ids)
        logging.info("Working with {} individuals".format(len(all_ids)))
        return all_ids

    def _load_pheno(self, individuals):
        """
        What tensorqtl wants: pandas DataFrame with pheno names as rows and
            individual IDs as columns
        """
        df = Parquet._read(self.pheno, to_pandas=True,
                           specific_individuals=individuals)
        logging.info("Loaded {} phenotypes".format(df.shape[1]))
        df['individual'] = df['individual'].astype('str')
        df = df.set_index('individual')
        self.covar_df = pandas.DataFrame(index=df.index)

        return df.T

def run(args):
    """
    1. Read plink and parquet files and find intersection of individuals.
    2. Load parquet phenotypes and get into correct data structure.
    3. Do some type of batching?
    4. run tensorqtl
    5. format results
    6. write results
    :param args:
    :return:
    """
    start = timer()
    file_out = FileOut(args.output_dir,
                       args.chromosome,
                       args.parquet_genotype_metadata)
    file_in = FileIn(args.plink_genotype,
                     args.parquet_phenotype,
                     args.rewrite_bim)
    pheno_df = file_in.get_pheno()
    covar_df = file_in.get_covar()
    geno_df = file_in.get_geno()
    logging.info("Computing summary statistics with p-val filter {} "
                 "and MAF filter {}".format(args.pval_filter, args.maf_filter))
    ss_df = trans.map_trans(geno_df,
                               pheno_df,
                               covar_df,
                               batch_size=2000,
                               pval_threshold=args.pval_filter,
                               maf_threshold=args.maf_filter,
                               verbose=False)
    gwas_df = file_out.make_gwas(ss_df, len(file_in.individuals))
    n_genes = len(gwas_df['gene'].drop_duplicates())
    for i, (pheno, df_i) in enumerate(gwas_df.groupby('gene')):
        logging.log(7, "Writing pheno {}/{}".format(i, n_genes))
        file_out.write_text(df_i, pheno)

    logging.info("Finished in {:.2f} seconds".format(timer() - start))






if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser("Fast GWAS summary statistics from tensorQTL")

    parser.add_argument("-plink_genotype")
    parser.add_argument("--parquet_genotype_metadata")
    parser.add_argument("--rewrite_bim", default=False, action='store_true')
    parser.add_argument("-parquet_phenotype")
    parser.add_argument("-output_dir")
    parser.add_argument("-chromosome")
    parser.add_argument('--parsimony', default=10, type=int)
    parser.add_argument('--pval_filter', type=float, default=1.0)
    parser.add_argument('--maf_filter', type=float, default=0.05)

    args = parser.parse_args()

    Logging.configure_logging(level=args.parsimony, target=sys.stdout)

    run(args)
