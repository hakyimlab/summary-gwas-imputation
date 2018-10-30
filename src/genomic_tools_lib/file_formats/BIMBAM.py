__author__ = "alvaro barbeira"

from .. import Utilities
from ..individual_data import Genotype

def load_bimbam_mean(file, sep=None):
    raise RuntimeError("Not implemented")

def write_bimbam_mean(genotype, path, to_minor=False):
    def to_lines(genotype):
        if to_minor:
            genotype = Genotype._to_minor_allele_frequency(genotype)
        variants = genotype.variants
        metadata = genotype.metadata
        for i,row in enumerate(metadata.itertuples()):
            minor_allele = row.allele_0
            major_allele = row.allele_1
            dosage = list(map(str,variants[i]))
            line = " ".join([str(row.id), minor_allele, major_allele]+dosage) + "\n"
            yield line
    Utilities.write_to_file(path, to_lines(genotype))

def write_bimbam_snp_annotation(genotype, path, to_minor=False):
    def to_lines(genotype):
        if to_minor:
            genotype = Genotype._to_minor_allele_frequency(genotype)
        metadata = genotype.metadata
        for i,row in enumerate(metadata.itertuples()):
            line = "{} {} {}\n".format(row.id, row.position, row.chromosome)
            yield line
    Utilities.write_to_file(path, to_lines(genotype))

def write_bimbam_pheno(phenotype, path):
    Utilities.save_dataframe(phenotype, path, header=False)

def save_study(study, prefix):
    Utilities.ensure_requisite_folders(prefix)

    genotype = study.get_genotype()
    pheno = study.get_phenotype()

    geno_ = prefix + ".geno.txt.gz"
    write_bimbam_mean(genotype, geno_)

    snp_ = prefix + ".snp.txt"
    write_bimbam_snp_annotation(genotype, snp_)

    pheno_ = prefix + ".pheno.txt"
    write_bimbam_pheno(pheno, pheno_)