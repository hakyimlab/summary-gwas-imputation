# -----------------------------------------------------------------------------
# SETUP ENVIRONMENT

suppressPackageStartupMessages(suppressWarnings(library(MatrixEQTL)))
suppressPackageStartupMessages(suppressWarnings(library(dplyr)))


# suppressPackageStartupMessages(suppressWarnings(library(argparse)))
# suppressPackageStartupMessages(suppressWarnings(library(arrow)))
# suppressPackageStartupMessages(suppressWarnings(library(rlang)))

log_stdout <- function(message, date=TRUE){
  if (date){
    message <- paste(Sys.time(), message)
  }
  print(message, stdout())
}

# -----------------------------------------------------------------------------
# LOAD AND PROCESS REGIONS FILE

load_regions <- function(fp, chr_){
  df <- read.table(fp, header=TRUE, sep="\t")
  df$chromosome <- as.integer(gsub('chr', '', df$chr ))
  df <- (df %>% filter(chromosome == chr_))
  df$index <- 0:(nrow(df)-1)
  df$fp_suffix <- paste0(paste(gsub('\\s', '', df$chr), df$start, df$stop,
                        sep='_'), '.txt')
  return(df)
}

# mutate_cond <- function(.data, condition, ..., envir = parent.frame()) {
#   # MANY THANKS to StackOverflow user G. Grothendieck for this code.
#   # https://stackoverflow.com/questions/34096162/dplyr-mutate-replace-several-columns-on-a-subset-of-rows
#   condition <- eval(substitute(condition), .data, envir)q
#   .data[condition, ] <- .data[condition, ] %>% mutate(...)
#   .data
# }
#
# annotations_to_regions <- function(regions, annotations){
#   # RETURNS a list with names of the index regions and values vectors of
#   #   applicaple variantIDs
#   blocks <- data.frame(annotations)
#   blocks$region_id <- rep_along(blocks$variant_id, NA)
#   for ( i in regions$index ){
#     start_ <- regions[i+1, 'start']
#     stop_ <- regions[i+1, 'stop']
#     blocks <- (blocks %>%
#                  # filter( pos >= start_ & pos <= stop_) %>%
#                  # mutate(region = i)
#                mutate_cond(position >= start_ & position <= stop_,
#                             region_id = i))
#   }
#   return(blocks)
# }

# -----------------------------------------------------------------------------
# GENOTYPE LOADING FUNCTION


# get_genotype <- function(fp, eids, cols=NULL){
#     if (is.null(cols)){
#       geno_df <- data.frame(read_parquet(fp))
#     }
#     else{
#       geno_df <- data.frame(read_parquet(fp, col_select=cols))
#     }
#     geno_df <- geno_df %>% filter(individual %in% eids)
#     rownames(geno_df) <- geno_df$individual
#     geno_df <- (geno_df %>%
#         select(-c(individual)) %>%
#         data.matrix %>%
#         t)
#     genos <- SlicedData$new()
#     genos$CreateFromMatrix(geno_df)
#     return(genos)
# }

matrix_to_sliceddata <- function(m, transpose=FALSE){
  if (transpose) {
  m <- t(m)
  }
  s_d <- SlicedData$new()
  s_d$CreateFromMatrix(m)
  return(s_d)
}

# -----------------------------------------------------------------------------
# LOAD SNP ANNOTATIONS
# GWAS_COLS <- c('variant_id', 'panel_variant_id', 'chromosome',
#                           'position', 'effect_allele', 'non_effect_allele',
#                           'current_build', 'frequency', 'sample_size', 'zscore',
#                           'pvalue', 'effect_size', 'standard_error',
#                           'imputation_status', 'n_cases', 'region_id')
#
# load_annotations <- function(fp){
#   snp_annotations <- read_parquet(fp)
#   na_cols <- c('imputation_status', 'n_cases')
#   for ( i in na_cols ){
#     snp_annotations[[i]] <- rep_along(snp_annotations$position, 'NA')
#   }
#   snp_annotations[['n_cases']] <- rep_along(snp_annotations$position, '10648')
#   snp_annotations[['current_build']] <- rep_along(snp_annotations$position,
#                                                   'hg19')
#   rename_cols <- c('non_effect_allele'='allele_0',
#                     'effect_allele'='allele_1',
#                     'frequency'='allele_1_frequency',
#                     'variant_id'='id')
#   snp_annotations <- snp_annotations %>% rename( rename_cols)
#   return(snp_annotations)
# }

# -----------------------------------------------------------------------------
# LOAD PHENOTYPE FROM PARQUET AND MAP COLUMN NUMBERS

# phenos <- <SlicedData object> Individuals on columns and phenos (genes) on
#           rows
# load_phenos <- function(fp, pheno_map){
#   phenos <- data.frame(read_parquet(fp))
#   pheno_map <- read.table(pheno_map, header=TRUE, sep="\t")
#   map_rearranged <- (pheno_map %>%
#                       slice(match(colnames(pheno_map), pheno_name)))
#   colnames(phenos) <- map_rearranged$pheno_num
#   rownames(phenos) <- phenos$individual
#   pheno_eids <- phenos$individual
#   phenos <- (phenos %>%
#       select(-c(individual)) %>%
#       data.matrix %>%
#       t)
#   pheno_data <- SlicedData$new()
#   pheno_data$CreateFromMatrix(phenos)
#   return(list(pheno_data, pheno_eids))
# }

# -----------------------------------------------------------------------------
# WRITE RESULTS BY IDP

# write_by_pheno <- function(results, out_dir, fp_suffix){
#   for (j in unique(results$gene)){
#     df_j <- results %>% filter(gene == j)
#     out_fp <- file.path(out_dir, paste('MatrixEQTL', j, fp_suffix, sep='_'))
#     write.table(df_j, file=out_fp, sep="\t", quote=FALSE, row.names=FALSE)
#   }
# }


# -----------------------------------------------------------------------------
# MAIN FUNCTIONS (CALLED BY PYTHON)
#
# init_data <- function(pheno, chr, geno, pheno_map){
#   init_lst <- list()
#   chr <- as.integer(chr)
#   init_lst[['chr']] <- chr
#   pheno_lst <- load_phenos(pheno)
#   init_lst[['pheno_data']] <- pheno_lst[[1]]
#   init_lst[['pheno_eids']] <- pheno_lst[[2]]
#   init_lst[['geno_fp']] <- geno
#   return(init_lst)
# }

run_summ_stats <- function(geno_matrix, pheno_matrix){
  # geno_matrix: matrix. Individuals in rows
  # pheno_matrix: matrix. Individuals in rows

  #print(paste("Rows: ", nrow(geno_matrix), nrow(pheno_matrix), sep = " "))
  #print(paste("Cols: ", ncol(geno_matrix), ncol(geno_matrix), sep = " "))
  genos <- matrix_to_sliceddata(geno_matrix, transpose=TRUE)
  phenos <- matrix_to_sliceddata(pheno_matrix, transpose=TRUE)
  me_new <- Matrix_eQTL_main(genos,
                               phenos,
                               output_file_name=NULL,
                               pvOutputThreshold=1.0,
                               useModel=modelLINEAR,
                               verbose=FALSE)
  return(me_new$all$eqtls)
}

# write_ <- function(df, path){
#   write.table(df, path, quote=FALSE, sep="\t", row.name=FALSE)
# }
#
# writer <- function(df, dir, chr){
#   phenos_ <- unique(df$gene) #TODO: Verify this
#   for (idp in idps_) {
#     df_i <- df %>% filter(gene == idp)
#     fname <- paste0('MatrixEQTL_', idp, '_chr', chr, '.txt')
#     fp <- file.path(dir, fname)
#     write_(df_i, fp)
#   }
# }
#
# GWAS_COLS <- c('variant_id', 'panel_variant_id', 'chromosome',
#                           'position', 'effect_allele', 'non_effect_allele',
#                           'current_build', 'frequency', 'sample_size', 'zscore',
#                           'pvalue', 'effect_size', 'standard_error',
#                           'imputation_status', 'n_cases', 'region_id')
#
# run <- function(args){
#   annotations <- annotations_to_regions(load_regions(args$regions, args$chr),
#                           load_annotations(args$annotations))
#   annotations <- annotations %>% filter(! is.na(region_id))
#   geno_pheno_data <- init_data(args$pheno, args$chr, args$geno, args$pheno_map)
#   summ_stats <- run_summ_stats(geno_pheno_data, annotations$variant_id)
#   summ_stats$zscore <- summ_stats$beta / summ_stats$statistic
#   rename_summ_stats <- c('variant_id'='snps', 'effect_size'='beta')
#   summ_stats <- summ_stats %>% rename(rename_summ_stats)
#   joined_df <- left_join(summ_stats, annotations, by=c('snp'='variant_id'))
#   writer(joined_df, args$out_dir, args$chr)
# }
#
# parser <- ArgumentParser()
# parser$add_argument('-chr')
# parser$add_argument('-regions')
# parser$add_argument('-annotations')
# parser$add_argument('-geno')
# parser$add_argument('-pheno')
# parser$add_argument('-pheno_map')
#
# args <- parser$parse_args()
#
# run(args)
