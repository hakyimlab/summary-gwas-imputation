__author__ = "alvaro barbeira"


def add_gwas_arguments_to_parser(parser):

    parser.add_argument("-separator", help="Character or string separating fields in input file. Defaults to any whitespace.", default=None)

    parser.add_argument("-skip_until_header", default=None,
                        help="Some files may be malformed and contain unespecified bytes in the beggining."
                             " Specify this option (string value) to identify a header up to which file contents should be skipped.")

    parser.add_argument("--handle_empty_columns", default=False, action="store_true",
                    help="Some files have empty columns, with values not even coded to missing. This instructs the parser to handle those lines."
                         "Be sur eyou want to use this.")

    parser.add_argument("--force_special_handling", default=False, action="store_true", help="Use custom text reading anyways.")

    parser.add_argument("-input_pvalue_fix",
                        help="If input GWAS pvalues are too small to handle, replace with these significance level. Use -0- to disable this behaviour and discard problematic snps.",
                        type=int, default=1e-50)

    parser.add_argument("-fill_from_snp_info", help="Add missing columns where appropriate from snp info", default=None)

    parser.add_argument("--enforce_numeric_columns", help="Force conversion to numeric type for any of  -beta, or, se, zscore, pvalue- ", action="store_true")

def discard_ambiguous(gwas, copy=True):
    if copy:
        gwas = gwas.copy()
    gwas["_key"] = gwas.effect_allele.str.lower() + gwas.non_effect_allele.str.lower()
    gwas = gwas.loc[~((gwas._key == "at") | (gwas._key == "ta") | (gwas._key == "cg") | (gwas._key == "gc"))]
    return  gwas.drop("_key", axis=1)

def get_chromosome_position_tree(d):
    r={}
    for t in d.itertuples():
        chromosome = t.chromosome
        if type(chromosome) == int:
            chromosome = "chr{}".format(chromosome)
        if not chromosome in r:
            r[chromosome] = set()
        c = r[chromosome]
        c.add(str(t.position))
    return r

#TODO: support for variable
def get_filter(chromosome_position_tree):
    def _filter(comps, tree):
        chrom = comps[2]
        if not chrom in tree:
            return True
        c = tree[chrom]

        if not comps[3] in c:
            return True

        return False

    return (lambda comps: _filter(comps, chromosome_position_tree))
