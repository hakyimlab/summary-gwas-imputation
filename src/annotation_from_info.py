import logging
import os
import gzip
import pandas

from genomic_tools_lib import Utilities, Logging

def load_dbsnp_mapping(path, variants_key):
    m = {}
    with gzip.open(path) as f:
        for i,l in enumerate(f):
            # if i>20000:
            #     break
            comps = l.decode().strip().split()

            obs_alleles = comps[9].split("/")
            if len(obs_alleles) < 2:
                continue

            chr = comps[1]
            start = comps[2]

            rsid = comps[4]
            strand = comps[6]
            ref_allele = comps[7]
            var_type = comps[11]

            alt_alleles_ = [x for x in obs_alleles if x != ref_allele]
            alt_alleles = set(alt_alleles_)

            if "-" == ref_allele or "-" in alt_alleles:
                pass
            else:
                start = str(int(start)+1)
            k = "{}_{}".format(chr, start)

            if not k in variants_key:
                continue

            for variant in variants_key[k]:
                variant_comps = variant.split("_")
                ref = variant_comps[2]
                alt = variant_comps[3]
                if ref == ref_allele and alt in alt_alleles:
                    m[variant] = rsid
                elif ref in alt_alleles and alt == ref_allele:
                    m[variant] = rsid
                elif ref_allele == "-":  # insertion
                    for alt_ in alt_alleles:
                        if alt_ == alt[1:]:
                            m[variant] = rsid
                            break
                elif "-" in alt:  # deletion or in-del
                    if len(ref) < len(alt):
                        kept = ref[:len(alt)]
                        deleted = ref[-len(alt):]
                        if deleted in alt_alleles:
                            m[variant] = rsid

    return m


import re
r_ = re.compile(".*chr(\d+).*")
def chr_key(x):
    k = os.path.split(x)[1]
    return int(r_.search(x).group(1))

def run(args):
    Utilities.ensure_requisite_folders(args.output)

    logging.info("Loading db snp file")
    #db_snp_mapping = load_dbsnp_mapping(args.dbsnp_file)

    logging.info("processing")
    files = sorted(args.info_files, key=chr_key)
    r = []
    variant_key = {}
    for p in files:
        with gzip.open(p) as f:
            logging.info("%s", p)
            for i,l in enumerate(f):
                if i == 0:
                    continue
                # if i > 20000:
                #     break
                comps = l.decode().strip().split()
                variant = comps[0]

                variant_comps = variant.split(":")
                chr = "chr"+variant_comps[0]
                pos = variant_comps[1]
                ref =  variant_comps[2]
                alt = variant_comps[3]
                if "CN" in ref or "CN" in alt:
                    continue
                freq = comps[3]

                variant_id = "{}_{}_{}_{}_b37".format(chr, pos, ref, alt)
                r.append((chr, pos, variant_id, ref, alt, freq))

                k = "{}_{}".format(chr, pos)
                if not k in variant_key:
                    variant_key[k] = []
                variant_key[k].append(variant_id)

    r = pandas.DataFrame(data=r, columns=["chromosome", "position", "id", "allele_0", "allele_1", "allele_1_frequency"])

    variant_key = {}
    for t in r.itertuples():
        k = "{}_{}".format(t.chromosome, t.position)
        if not k in variant_key:
            variant_key[k] = []
        variant_key[k].append(t.id)

    logging.info("looking for rsids in ucsc dbsnp file")
    dbsnp_mapping = load_dbsnp_mapping(args.dbsnp_file, variant_key)

    rsids = []
    for id in r.id:
        rsids.append(dbsnp_mapping[id] if id in dbsnp_mapping else "NA")
    r["rsid"] = rsids
    logging.info("Saving")
    Utilities.save_dataframe(r, args.output)

    logging.info("Done")

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser("Convert -imputation info- file to variant metadata format")
    parser.add_argument("-info_files", nargs="+")
    parser.add_argument("-dbsnp_file")
    parser.add_argument("-parsimony", type=int, default=logging.INFO)
    parser.add_argument("-output")
    args = parser.parse_args()

    Logging.configure_logging(args.parsimony)

    run(args)