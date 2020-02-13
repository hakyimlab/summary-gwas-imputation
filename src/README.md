# summary-gwas-imputation/src

## Data Formats

### Model Training Format

Possibly g-zipped text file with tabs separating columns and a variant on each row, and an individual on each column:
```
varID   ID1 ID2 ID3 ... 
chr1_13526_C_T_b38  0   1   0   ...
```

## Directory

### model_training_genotype_to_parquet.py

Takes a genotype file in the [model training format](#Model-Training-Format) and creates a parquet binary file. 
