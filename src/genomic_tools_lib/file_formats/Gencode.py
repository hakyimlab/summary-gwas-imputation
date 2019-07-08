__author__ = "alvaro barbeira"
from .. import Utilities

K_NOT_GENES = {"transcript","exon","CDS","UTR","start_codon","stop_codon","Selenocysteine"};

# look at gencode http://www.gencodegenes.org/data_format.html
class GFTF:
    """gencode file table format"""
    CHROMOSOME = 0
    ANNOTATION_SOURCE = 1
    FEATURE_TYPE = 2
    N_START_LOCATION = 3
    N_END_LOCATION = 4
    SCORE = 5
    GENOMIC_STRAND = 6
    GENOMIC_PHASE = 7
    KEY_VALUE_PAIRS = 8

    #there are several other key-value pairs but we are concerned with these
    K_GENE_ID = "gene_id"
    K_TRANSCRIPT_ID = "transcript_id"
    K_GENE_TYPE = "gene_type"
    K_GENE_STATUS = "gene_status"
    K_GENE_NAME = "gene_name"
    K_TRANSCRIPT_TYPE = "transcript_type"
    K_TRANSCRIPT_STATUS = "transcript_status"
    K_TRANSCRIPT_NAME = "transcript_name"
    K_EXON_NUMBER = "exon_number"
    K_EXON_ID = "exon_id"
    K_LEVEL = "level"

    #some are missing
    TAG = "tag"

def load(path, gene_ids=None, feature_type_whitelist={"gene"}, gene_type_white_list=None, transcript_type_whitelist=None, selected_key_value_pairs=[GFTF.K_GENE_ID, GFTF.K_GENE_NAME], collapse_strand=False):
    results = []

    F = GFTF
    t_ = str.maketrans(dict.fromkeys(';"'))
    for i,line in Utilities.iterate_file(path):
        if "#" == line[0]: continue

        comps = tuple(line.strip().split())
        feature_type = comps[F.FEATURE_TYPE]

        if feature_type_whitelist and feature_type not in feature_type_whitelist: continue

        key_value_pairs = [x.translate(t_) for x in comps[GFTF.KEY_VALUE_PAIRS:]]
        key_value_pairs = {key_value_pairs[i]:key_value_pairs[i+1] for i in range(0, len(key_value_pairs), 2)}

        if gene_ids:
            gene_id = key_value_pairs[GFTF.K_GENE_ID]
            if not gene_id in gene_ids: continue

        if transcript_type_whitelist:
            if not GFTF.K_TRANSCRIPT_TYPE in key_value_pairs: continue
            transcript_type = key_value_pairs[GFTF.K_TRANSCRIPT_TYPE]
            if not transcript_type in transcript_type_whitelist: continue

        if gene_type_white_list:
            gene_type = key_value_pairs[GFTF.K_GENE_TYPE]
            if not gene_type in gene_type_white_list: continue

        if not collapse_strand:
            r = (comps[F.CHROMOSOME], comps[F.N_START_LOCATION], comps[F.N_END_LOCATION], feature_type, comps[F.GENOMIC_STRAND])
        else:
            if comps[F.GENOMIC_STRAND] == "+":
                r = (comps[F.CHROMOSOME], comps[F.N_START_LOCATION], comps[F.N_END_LOCATION], feature_type, comps[F.GENOMIC_STRAND])
            else:
                r = (comps[F.CHROMOSOME], comps[F.N_END_LOCATION], comps[F.N_START_LOCATION], feature_type, comps[F.GENOMIC_STRAND])

        r += tuple([key_value_pairs[x] for x in selected_key_value_pairs])
        results.append(r)

    columns = ["chromosome", "start_location", "end_location", "feature_type", "strand"]+selected_key_value_pairs
    results = Utilities.to_dataframe(results, columns)
    return results
