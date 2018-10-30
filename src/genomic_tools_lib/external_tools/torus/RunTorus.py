__author__ = "alvaro barbeira"
import os
import re
import gzip
import shutil
import logging
import traceback

from subprocess import call

import pyarrow.parquet as pq

from ... import Exceptions
from ... import Utilities
from ...data_management import TextFileTools
from ...file_formats import Gencode
from ...file_formats.eqtl import GTEx
from ...Exceptions import ReportableException

class _Context():
    def get_input_eqtl(self): raise Exceptions.ReportableException("Not implemented")
    def get_variant_whitelist(self): raise Exceptions.ReportableException("Not implemented")
    def get_input_gene_annotation(self): raise Exceptions.ReportableException("Not implemented")
    def get_intermediate_folder(self): raise Exceptions.ReportableException("Not implemented")
    def get_torus_exe(self): raise Exceptions.ReportableException("Not implemented")
    def get_output_folder(self): raise Exceptions.ReportableException("Not implemented")
    def get_delete_intermediate(self): raise Exceptions.ReportableException("Not implemented")

def _pq_from_eqtl_to_torus(input_path, output_path, variant_whitelist):
    snps=set()
    genes=set()
    with TextFileTools.TextDataFrameSink(output_path, write_header=False) as sink:
        p = pq.ParquetFile(input_path)
        for i in range(0, p.num_row_groups):
            logging.log(8, "Processing gene %d/%d", i+1, p.num_row_groups)
            d = p.read_row_group(i).to_pandas()
            d = d.assign(zscore = d.beta/d.se)
            d = d[["variant_id", "gene", "zscore"]]

            if variant_whitelist:
                #help pandas a bit with the isin() method
                w = {x for x in d.variant_id if x in variant_whitelist}
                d = d.loc[d.variant_id.isin(w)]

            if d.shape[0] == 0:
                continue

            sink.sink(d)
            genes.add(d.gene.values[0])
            snps.update(d.variant_id.values)
    return snps, genes

def _t_from_eqtl_to_torus(input_path, output_path, variant_whitelist):
    snps=set()
    genes=set()
    with gzip.open(output_path, "w") as o:
        g, v, b, s = 0, 1, 7, 8
        for i,line in Utilities.iterate_file(input_path, skip_first=True):
            comps = line.strip().split()
            _v = comps[v]
            if not _v in variant_whitelist:
                continue
            _g = comps[g]
            _z = float(comps[b])/float(comps[s])
            line = "{} {} {:.2f}\n".format(_v, _g, _z)
            o.write(line.encode())
            snps.add(_v)
            genes.add(_g)
    return  snps, genes


def from_eqtl_to_torus(input_path, output_path, variant_whitelist):
    if re.search(".parquet$", input_path):
        return _pq_from_eqtl_to_torus(input_path, output_path, variant_whitelist)
    elif re.search(".allpairs.txt.gz$", input_path):
        return _t_from_eqtl_to_torus(input_path, output_path, variant_whitelist)
    else:
        raise  RuntimeError("Unsupported file")

def generate_torus_snp_map(snps, snp_path):
    with gzip.open(snp_path, "w") as f:
        f.write("SNP\tChromosome\tPosition\n".encode())
        for r in snps:
            comps = r.split("_")
            chr = comps[0]
            pos = comps[1]
            if "X" in chr:
                continue
            chr = chr.split("chr")[1]
            f.write( Utilities.to_line([r, chr, pos]).encode() )

def from_gene_annotation_to_torus(input_path, genes, output_path):
    annotation = Gencode.load(input_path, gene_ids=genes)
    annotation = annotation.rename(columns={"chromosome":"Chromosome", "gene_id":"Gene", "start_location":"TSS", "end_location":"TES"})
    annotation = annotation[["Gene","Chromosome", "TSS", "TES"]]
    annotation.Chromosome = annotation.Chromosome.str.split("chr").str.get(1)
    Utilities.save_dataframe(annotation, output_path, header=False)

def run_torus(context):
    try:
        _run_torus(context)
    except ReportableException as ex:
        logging.info("Reportable exception running torus: %s", ex.msg)
    except Exception as ex:
        logging.info("Exception running torus:\n%s", traceback.format_exc())
    finally:
        if context.get_delete_intermediate():
            folder = context.get_intermediate_folder()
            if os.path.exists(folder):
                shutil.rmtree(folder)

def _run_torus(context):
    Utilities.maybe_create_folder(context.get_intermediate_folder())
    Utilities.maybe_create_folder(context.get_output_folder())

    intermediate = context.get_intermediate_folder()

    logging.info("Processing eqtl file")
    eqtl_path = os.path.join(intermediate, "eqtl.gz")

    snps, genes = from_eqtl_to_torus(context.get_input_eqtl(), eqtl_path, context.get_variant_whitelist())

    logging.info("Processing snp map")
    snp_map_path = os.path.join(intermediate, "snp.gz")
    generate_torus_snp_map(snps, snp_map_path)

    logging.info("Processing gene annotation")
    gene_annotation_path = os.path.join(intermediate, "gene.gz")
    from_gene_annotation_to_torus(context.get_input_gene_annotation(), genes, gene_annotation_path)

    output_priors = os.path.join(context.get_output_folder(), "priors")

    command = _torus_command(context.get_torus_exe(), eqtl_path, snp_map_path, gene_annotation_path, output_priors).split()

    output_estimate = os.path.join(context.get_output_folder(), "estimate.txt")
    output_log = os.path.join(context.get_output_folder(), "torus_log.txt")

    with open(output_estimate, "w") as o:
        with open(output_log, "w") as e:
            logging.info("Runing torus")
            logging.log(9, "command: '%s'", command)
            call(command, stdout=o, stderr=e)

def _torus_command(torus_exe, eqtl, snp, gene, prior, _call=True):
    command=\
r"""{} \
-d {} --load_zval \
-smap {} \
-gmap {} \
-dump_prior {}
""".format(torus_exe, eqtl, snp, gene, prior)
    if _call:
        command=command.replace("\\", "")
    return command