import logging
import os
import gzip
import re
from timeit import default_timer as timer

import numpy
import pandas
import pyliftover

from genomic_tools_lib import Logging, Utilities
from genomic_tools_lib.file_formats import DBSnp
from genomic_tools_lib.file_formats.gwas import GWAS, Utilities as GWASUtilities

import gwas_parsing

########################################################################################################################
class alelle_check:
    def matches(self, panel_ref_allele, panel_alt_allele, ref_allele, alt_alleles):
        raise RuntimeError("Not implemented")

    def parse(self, panel_ref_allele, panel_alt_allele, ref_allele, alt_alleles):
        raise RuntimeError("Not implemented")

class snp_strand_p1_swap_p1(alelle_check):
    def matches(self, panel_ref_allele, panel_alt_allele, ref_allele, alt_alleles):
        return panel_ref_allele == ref_allele and panel_alt_allele in alt_alleles

    def parse(self, panel_ref_allele, panel_alt_allele, ref_allele, alt_alleles):
        return  1, 1, panel_ref_allele, panel_alt_allele

class snp_strand_p1_swap_m1(alelle_check):
    def matches(self, panel_ref_allele, panel_alt_allele, ref_allele, alt_alleles):
        return panel_ref_allele in alt_allele and panel_alt_allele == ref_alleles

    def parse(self, panel_ref_allele, panel_alt_allele, ref_allele, alt_alleles):
        return  1, -1, panel_alt_allele, panel_ref_allele


########################################################################################################################

def l_(comps):
    return "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(*comps).encode()


def run(args):
    if os.path.exists(args.output):
        logging.info("Output already exists, nope.")
        return

    Utilities.ensure_requisite_folders(args.output)
    Utilities.ensure_requisite_folders(args.discard)

    if args.liftover:
        logging.info("Acquiring liftover")
        l = pyliftover.LiftOver(args.liftover)
    else:
        logging.info("Will not perform lift over")
        l = None

    logging.info("Loading snp reference metadata")
    snp_reference_metadata = pandas.read_table(args.snp_reference_metadata)
    reference = {}
    for t in snp_reference_metadata.itertuples():
        k = "chr{}_{}".format(t.chromosome, t.position)
        if k in reference:
            raise RuntimeError("coordinate is already present")
        reference[k] = (t.id, t.rsid)

    dbsnp_format = {x: i for i, x in enumerate(DBSnp.DBSNP._fields)}
    complement_translation = "CGTA".maketrans({"C": "G", "G": "C", "T":"A", "A": "T"})

    logging.info("Processing db snp file")
    if args.discard:
        discard = gzip.open(args.discard, "w")
        discard.write(l_(["rsid", "chromosome", "position", "a0", "a1", "strand", "type", "panel_variant_id", "panel_variant_rsid", "panel_variant_a0", "panel_variant_a1", "swap", "strand_reversal"]))

    allele_re = re.compile("chr\d+_\d+_(.*)_(.*)_b38")

    with gzip.open(args.output, "w") as result:
        result.write(l_(["rsid", "chromosome", "position", "a0", "a1", "strand", "type", "panel_variant_id", "panel_variant_rsid", "panel_variant_a0", "panel_variant_a1", "swap", "strand_reversal"]))
        with gzip.open(args.db_snp_file) as db_snp:
            db_snp.readline()
            for i,line in enumerate(db_snp):
                comps = line.decode().strip().split("\t")

                obs_alleles = comps[9].split("/")
                if len(obs_alleles) < 2:
                    continue

                chr = comps[1]
                start_0 = comps[2]
                _new_chromosome, _new_position = gwas_parsing._lift(l, chr, start_0) if l else (chr, int(start_0))

                if _new_chromosome == "NA" or _new_position == "NA":
                    continue

                k = "{}_{}".format(_new_chromosome, _new_position+1)
                if not k in reference:
                    continue

                rsid = comps[4]
                strand = comps[6]
                ref_allele = comps[7]
                var_type = comps[11]

                alt_alleles_ = [x for x in obs_alleles if x != ref_allele]
                alt_alleles = set(alt_alleles_)

                panel_variant_id, panel_variant_rsid = reference[k]
                panel_variant_rsid = panel_variant_rsid if type(panel_variant_rsid) == str else "NA"
                panel_alleles = allele_re.search(panel_variant_id)
                panel_ref_allele = panel_alleles.group(1)
                panel_alt_allele = panel_alleles.group(2)

                strand_reversed_panel_ref_allele = panel_ref_allele.translate(complement_translation)
                strand_reversed_panel_alt_allele = panel_alt_allele.translate(complement_translation)
                # if args.reverse_swap:
                #     strand_reversed_panel_ref_allele = strand_reversed_panel_ref_allele[::-1]
                #     strand_reversed_panel_alt_allele = strand_reversed_panel_alt_allele[::-1]

                swap, strand_reversal, selected_ref_allele, selected_alt_allele = None, None, ref_allele, alt_alleles_[0]
                if len(panel_ref_allele) == 1 and len(panel_alt_allele) == 1:
                    #snp
                    if panel_ref_allele == ref_allele and panel_alt_allele in alt_alleles:
                        swap, strand_reversal, selected_ref_allele, selected_alt_allele =  1,  1, panel_ref_allele, panel_alt_allele
                    elif panel_ref_allele in alt_alleles and panel_alt_allele == ref_allele:
                        swap, strand_reversal, selected_ref_allele, selected_alt_allele = -1,  1, panel_alt_allele, panel_ref_allele
                    elif strand_reversed_panel_ref_allele == ref_allele and strand_reversed_panel_alt_allele in alt_alleles:
                        swap, strand_reversal, selected_ref_allele, selected_alt_allele =  1, -1, strand_reversed_panel_ref_allele, strand_reversed_panel_alt_allele
                    elif strand_reversed_panel_ref_allele in alt_alleles and strand_reversed_panel_alt_allele == ref_allele:
                        swap, strand_reversal, selected_ref_allele, selected_alt_allele = -1, -1, strand_reversed_panel_alt_allele, strand_reversed_panel_ref_allele
                elif len(panel_ref_allele) > 1 and len(panel_alt_allele) == 1 and ref_allele != "-":
                    #deletion
                    deleted = panel_ref_allele[1:]
                    strand_reversed_deleted = strand_reversed_panel_ref_allele[1:]
                    # if args.reverse_swap:
                    #     strand_reversed_deleted = strand_reversed_panel_ref_allele[:-1]
                    for si_, allele_ in enumerate(alt_alleles):
                        if allele_ == deleted:
                            swap, strand_reversal, selected_ref_allele, selected_alt_allele =  1,  1, allele_, "-"
                        if allele_ == strand_reversed_deleted:
                            swap, strand_reversal, selected_ref_allele, selected_alt_allele =  1, -1, allele_, "-"
                elif len(panel_ref_allele) == 1 and len(panel_alt_allele) > 1 and ref_allele == "-":
                    inserted = panel_alt_allele[1:]
                    strand_reversed_inserted = strand_reversed_panel_alt_allele[1:]#[:-1]
                    # if args.reverse_swap:
                    #     strand_reversed_inserted = strand_reversed_panel_alt_allele[:-1]
                    for si_, allele_ in enumerate(alt_alleles):
                        if allele_ == inserted:
                            swap, strand_reversal, selected_ref_allele, selected_alt_allele =  1,  1, "-", allele_
                        if allele_ == strand_reversed_inserted:
                            swap, strand_reversal, selected_ref_allele, selected_alt_allele =  1, -1, "-", allele_

                else:
                    pass

                ol = l_([rsid, chr, str(int(start_0) + 1), selected_ref_allele, selected_alt_allele, strand, var_type, panel_variant_id, panel_variant_rsid, panel_ref_allele, panel_alt_allele, swap, strand_reversal])
                if swap is not None and strand is not None and selected_ref_allele is not None and selected_alt_allele is not None:
                    result.write(ol)
                else:
                    discard.write(ol)
    discard.close()
    logging.info("Done")

if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser("parse a GWAS into a standard format")
    parser.add_argument("-db_snp_file", help="snp definition from a particular hg release")
    parser.add_argument("-liftover", help = "File with liftover chain")
    parser.add_argument("-snp_reference_metadata", help="File with reference (GTEx) snp metadata")
    parser.add_argument("-output", help="Where the output should go")
    parser.add_argument("-discard", help="Where the discard output should go")
    parser.add_argument("-verbosity", help="Log verbosity level. 1 is everything being logged. 10 is only high level messages, above 10 will hardly log anything", default = 10, type = int)
    parser.add_argument("--reverse_swap", action="store_true")
    GWASUtilities.add_gwas_arguments_to_parser(parser)
    args = parser.parse_args()

    Logging.configure_logging(args.verbosity)

    run(args)