# warning

`src/compute_genomic_mapping.py` is a heuristic tool to map variants from an older reference to a new one.
It is a honest, if futile, attempt at covering a frequent or semi-frequent variant definitions
from a wild zoology of edge cases in UCSC variant definitions,
and map them to GTEx variant definition.

It doesn't handle all possible cases but should handle most.

## Edge case 1
See the following GTEx variant:

```
chr   variant_pos variant_id                ref alt  num_alt_per_site	rs_id_dbSNP150_GRCh38p7
...
chr10 47225796    chr10_47225796_G_GTCC_b38 G   GTCC 1                  rs142868522

```  

This variant explictly declares alleles `G` and `GTCC` (rs142868522 in b38) 
but dbSNP declares alelles `G` and `GTCCTC`.
IN snp125_hg17, there is a variant that, after liftover, maps to GTEx definition: `rs10631152`
but it has been delisted. What to do then, in this cases?

## Edge case 2

Variant `rs11348115` in UCSC https://www.ncbi.nlm.nih.gov/snp/rs11348115
is listed as `-/T` (delT) but is displayes as deletion of `TTTTTTTTTTTTTTTT`

# Edge case 3

Many snps have multiple entries in the UCSC file: see for example rs3007516.

# Edge case 4

Both `rs9794226` and `rs11813951` have the same position, and have the same alleles after swapping strands.
