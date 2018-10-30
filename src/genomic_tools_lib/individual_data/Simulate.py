__author__ = "alvaro barbeira"
import numpy
import pandas
from . import Genotype, Study

_DEFAULT_SAMPLE_SIZE=500
_DEFAULT_VARIANTS_PER_CHROM=100

def reset_seed():
    numpy.random.seed(0)

def random_alleles():
    _map = {
        "A": ["C","G"],
        "T": ["C","G"],
        "C": ["A", "T"],
        "G": ["A", "T"]}
    first = numpy.random.choice(list(_map.keys()))
    second = numpy.random.choice(_map[first])
    return first, second

def random_dosage(sample_size):
    a = list(numpy.random.uniform(0,2.0/3, sample_size))
    b = list(numpy.random.uniform(2.0/3, 4.0/3, sample_size))
    c = list(numpy.random.uniform(4.0/3, 2, sample_size))
    mode_ = numpy.random.choice(3)
    if mode_ == 0:
        pool = a*5 + b + c
    elif mode_ == 1:
        pool = a + b*5 + c
    else:
        pool = a + b + c*5
    return numpy.random.choice(pool, sample_size, False)

def simulate_genotype(variants_per_chromosome=_DEFAULT_VARIANTS_PER_CHROM, sample_size=_DEFAULT_SAMPLE_SIZE):
    id_ = 0

    dosages = []
    metadata = []
    for i in range(1,23):
        chromosome = i
        pos_ = 1
        for j in range(0, variants_per_chromosome):
            id_ += 1
            rsid = "rs{}".format(id_)
            pos_ += numpy.random.randint(0,100)
            allele_0, allele_1 = random_alleles()
            dosage = random_dosage(sample_size)
            frequency = numpy.mean(dosage)/2

            metadata.append((chromosome, pos_, rsid, allele_0, allele_1, frequency))
            dosages.append(dosage)

    metadata = Genotype._metadata_from_raw(metadata)

    return Genotype.Genotype(dosages, metadata)

def simulate_individuals(sample_size=_DEFAULT_SAMPLE_SIZE):
    return ["ID_{}".format(i) for i in range(0,sample_size)]

def simulate_random_phenotype(individual_ids):
    p = numpy.random.uniform(size=len(individual_ids))
    return pandas.DataFrame({"phenotype": p})

def simulate_bslmm_phenotype(genotype, individual_ids, selected):
    sample_size = len(individual_ids)
    m = genotype.metadata
    d = genotype.variants
    pheno_b = numpy.zeros(sample_size)
    for i,row in enumerate(m.itertuples()):
        if row.id in selected:
            coef = numpy.random.uniform(1,2, sample_size)*numpy.random.choice([1,-1])
        else:
            coef = numpy.random.uniform(-0.02, 0.02, sample_size)
        d_ = numpy.nan_to_num(d[i])
        pheno_b += d_*coef
    pheno_b += numpy.random.normal(size=sample_size)
    pheno_r = numpy.random.normal(size=sample_size)
    pheno = pandas.DataFrame({"GID1":pheno_b, "GID2":pheno_r})
    return pheno

def simulate_bslmm_study(snps_per_chromosome=_DEFAULT_VARIANTS_PER_CHROM):
    snps_per_chromosome = snps_per_chromosome if snps_per_chromosome is not None else _DEFAULT_VARIANTS_PER_CHROM
    genotype = simulate_genotype(snps_per_chromosome)
    individuals = simulate_individuals()
    gene_annotation = simulate_gencode_parsed()
    selected_snps = select_snps(gene_annotation,  genotype.metadata)
    phenotype = simulate_bslmm_phenotype(genotype, individuals, selected_snps)
    return Study.Study(genotype, phenotype, individuals), selected_snps, gene_annotation

def simulate_gencode_parsed():
    d = [[1, "GID1", "A", 100, 300, "pseudogene"],
        [2, "GID2", "B", 4000, 5000, "lincRNA"]]
    d = pandas.DataFrame(d, columns=["chr", "gene_id", "gene_name", "start", "end", "gene_type"])
    return d

def select_snps(gene_annotation, metadata):
    g_ = gene_annotation[gene_annotation.gene_id == "GID1"]
    c_ = g_.chr[0]
    s_ = 0
    e_ = g_.end[0] + 100
    v_ = metadata[(metadata.chromosome == c_) & (metadata.position >= s_) & (metadata.position <= e_)]
    selected_ = numpy.random.choice(range(0, v_.shape[0]), 2)
    v_ = v_.iloc[selected_]
    return set(v_.id)