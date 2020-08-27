__author__ = "alvaro barbeira"
import sys
import gzip
import re
import numpy
import logging
import argparse
import os
from genomic_tools_lib import Logging

def _fix_int(dosage, flip):
    m = float(numpy.mean([int(x) for x in dosage if x != "NA"]))
    dosage = [int(x) if x != "NA" else m for x in dosage]
    if flip:
        dosage = map(lambda x: str(2-(x)), dosage)
    return list(map(str,dosage))

def _fix_dosage(dosage, flip):
    m = float(numpy.mean([float(x) for x in dosage if x != "NA"]))
    dosage = [float(x) if x != "NA" else m for x in dosage]
    if flip:
        dosage = map(lambda x: str(2-(x)), dosage)
    return list(map(str,dosage))

CHR=0
SNP=1
_C_M=2
POS=3
COUNTED=4
ALT=5
FIRST=6

def get_key(id_mapping):
    logging.info("Reading id file from {}".format(id_mapping))
    key = {}
    with gzip.open(id_mapping) as r:
        r.readline()
        for l in r:
            comps = l.strip().split()
            # from IPython import embed; embed(); exit()
            F = None
            if comps[2] == "TRUE":
                F = True
            if comps[2] == "FALSE":
                F = False
            key[comps[0]] = (comps[1], F)
    return key

def _write_with_key(i, line, a, d, did, count_f, key=None,
                    whitelist=None):
    c = line.strip().split()
    chr = c[CHR]


    var = c[SNP]
    if not var in key:
        return

    _k = key[var]
    if _k[1] is None:
        return

    if whitelist and not _k[0] in whitelist:
        return

    ref = c[ALT] if not _k[1] else c[COUNTED]
    alt = c[COUNTED] if not _k[1] else c[ALT]
    dosage = count_f(c[FIRST:], _k[1])
    rsid = _k[0]

    # allele names are in different convention
    var_id = "_".join(["chr" + chr, c[POS], ref, alt])
    if var_id in did:
        return
    did.add(var_id)

    a_ = "\t".join(["chr" + chr, c[POS], var_id, ref, alt, "NA", "NA", "NA", rsid]) + "\n"
    a.write(a_)

    d_ = "\t".join([var_id] + dosage) + "\n"
    d.write(d_)

def _write_without_key(i, line, a, d, did, count_f, key=None,
                    whitelist=None):
    c = line.strip().split()
    pos = c[POS]
    chr = c[CHR]
    rsid = c[SNP]
    ref = c[ALT]
    alt = c[COUNTED]
    dosage = count_f(c[FIRST:], False)
    var_id = "_".join(['chr'+ chr, pos, ref, alt])
    a_ = "\t".join(["chr" + chr, pos, var_id, ref, alt, "NA", "NA", "NA", rsid]) + "\n"
    a.write(a_.encode())
    d_ = "\t".join([var_id] + dosage) + "\n"
    d.write(d_.encode())

def run(args):

    if args.id_mapping:
        key = get_key(args.id_mapping)
        write_f = _write_with_key
    else:
        key = None
        write_f = _write_without_key

    if args.dosage:
        count_f = _fix_dosage
    else:
        count_f = _fix_int

    if args.whitelist:
        logging.info("Reading whitelist from {}".format(args.whitelist))
        mod="F_"
        whitelist=set()
        with open(args.whitelist) as w:
            for l in w:
                whitelist.add(l.strip())
    else:
        whitelist=None
        mod=""

    a_fp = args.output_prefix + '.a.txt.gz'
    if os.path.isfile(a_fp):
        raise ValueError("Annotation filepath already exists.")
    a = gzip.open(a_fp, 'w')

    d_fp = args.output_prefix + '.d.txt.gz'
    d = gzip.open(d_fp, 'w')


    did = set()
    logging.info("Converting {}".format(args.raw_plink))
    with open(args.raw_plink) as r:
        head = r.readline().strip().split()
        ids = [x.split("_")[1] for x in head[FIRST:]]
        header = "\t".join(["varID"] + ids) + "\n"
        a.write("chromosome\tpos\tvarID\tref_vcf\talt_vcf\tR2\tMAF\trsid\trsid_dbSNP150\n".encode())
        d.write(header.encode())
        for i, l in enumerate(r):
            write_f(i, l, a, d, did, count_f, key=key, whitelist=whitelist)
            if i % 100000 == 0:
                logging.debug('Finished with variant {}'.format(i))

    a.close()
    d.close()
    logging.info("Finished")

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-raw_plink')
    parser.add_argument('-output_prefix')
    parser.add_argument('--dosage', default=False, action='store_true')
    parser.add_argument('--id_mapping')
    parser.add_argument('--whitelist')
    parser.add_argument('--save_bad_keys')
    parser.add_argument('-parsimony', type=int, default=10)

    args = parser.parse_args()
    f = "%(asctime)s : %(levelname)s - %(message)s"

    Logging.configure_logging(level=args.parsimony, with_date=True, target=sys.stdout)

    run(args)
