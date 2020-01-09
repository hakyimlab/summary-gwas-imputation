import gzip
import os
import logging
from genomic_tools_lib import Utilities, Logging
from genomic_tools_lib.data_management import KeyedDataSource, TextFileTools

def _to_al(comps):
    return "{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(*comps).encode()

def _to_gl(comps, format):
    return format.format(*comps).encode()

def run(args):
    Utilities.ensure_requisite_folders(args.output_prefix)

    logging.info("Loading snp reference")
    key = KeyedDataSource.load_data(args.snp_reference_file, "variant_id", "rs_id_dbSNP150_GRCh38p7", value_conversion=KeyedDataSource.dot_to_na)
    logging.info("Loading samples")
    samples = TextFileTools.load_list(args.samples)
    genotype_format_string = "\t".join(["{}"]*(len(samples)+1))+ "\n"

    og = args.output_prefix + "_genotype.txt.gz"
    oa = args.output_prefix + "_annotation.txt.gz"
    if os.path.exists(og) or os.path.exists(oa):
        logging.info("Output exists. Nope.")
        return

    logging.info("Processing")
    with gzip.open(args.genotype) as geno:
        with gzip.open(og, "w") as _og:
            _og.write(_to_gl(["varID"]+samples, genotype_format_string))
            with gzip.open(oa, "w") as _oa:
                _oa.write(_to_al(["chromosome", "position", "id", "allele_0", "allele_1", "allele_1_frequency" , "rsid"]))
                for i, line in enumerate(geno):
                    comps = line.decode().strip().split()

                    chr = "chr"+ comps[0]
                    pos = comps[2]
                    ref = comps[3]
                    alt = comps[4]
                    af = comps[5]
                    dosage = comps[6:]

                    var_id = "{}_{}_{}_{}_b38".format(chr, pos, ref, alt)
                    if var_id in key:
                        id = key[var_id]
                        comps[1] = var_id
                        _og.write(_to_gl([var_id]+dosage, genotype_format_string))
                        _oa.write(_to_al([chr, pos, var_id, ref, alt, af, id]))
                        next

                    var_id = "{}_{}_{}_{}_b38".format(chr, pos, alt, ref)
                    if var_id in key and len(ref) == 1 and len(alt) == 1:
                        id = key[var_id]
                        af = str(1-float(af))
                        dosage = list(map(lambda x: str(2-int(x)), comps[6:]))
                        _og.write(_to_gl([var_id]+dosage, genotype_format_string))
                        _oa.write(_to_al([chr, pos, var_id, alt, ref, af, id]))
                        next

    logging.info("Finished conversion")

if __name__ == "__main__":
    import  argparse
    parser = argparse.ArgumentParser("Convert genotypes in -PrediXcan Format- to -Model training format-")
    parser.add_argument("-genotype")
    parser.add_argument("-samples")
    parser.add_argument("-snp_reference_file")
    parser.add_argument("-output_prefix")
    parser.add_argument("-parsimony", help="Log parsimony level. 1 is everything being logged. 10 is only high level messages, above 10 will hardly log anything", default=10, type =int)
    args = parser.parse_args()
    Logging.configure_logging(args.parsimony)

    run(args)
