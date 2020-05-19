suppressPackageStartupMessages(suppressWarnings(library(dplyr)))
suppressPackageStartupMessages(suppressWarnings(library(ggplot2)))

# -----------------------------------------------------------------------------
# LOGGING FUNCTION

log_stdout <- function(message, date=TRUE){
  if (date){
    message <- paste(Sys.time(), message)
  }
  print(message, stdout())
}

# -----------------------------------------------------------------------------
# MAKE PLOTS

#log_stdout("Adding functions to R namespace")
make_pval_plot <- function(pvalues, fp, pheno_name){
 # log_stdout("Opened make_pval_plot function")
  title_s <- paste0("Metabolite ", pheno_name)
  plt <- (gg_qqplot(pvalues)
            + labs(title = title_s)
            + scatter_base_theme_())
  #log_stdout("Beginning to write plot to the disk")
  png(file=fp)
  print(plt)
  dev.off()
  # ggsave(fp, height=9.76, width=9.76)

}



gg_qqplot <- function(ps, ci = 0.95) {
  # Many thanks to github user slowkow for this function
  # https://gist.github.com/slowkow/9041570
 # log_stdout("Opened gg_qqplot function")
  n  <- length(ps)
  df <- data.frame(
    observed = -log10(sort(ps)),
    expected = -log10(ppoints(n))
#    clower   = -log10(qbeta(p = (1 - ci) / 2, shape1 = 1:n, shape2 = n:1)),
#    cupper   = -log10(qbeta(p = (1 + ci) / 2, shape1 = 1:n, shape2 = n:1))
  )
  log10Pe <- expression(paste("Expected -log"[10], plain(P)))
  log10Po <- expression(paste("Observed -log"[10], plain(P)))
  ggplot(df) +
    geom_point(aes(expected, observed), size = 3) +
    geom_abline(intercept = 0, slope = 1, alpha = 0.5) +
#    geom_line(aes(expected, cupper), linetype = 2) +
#    geom_line(aes(expected, clower), linetype = 2) +
    xlab(log10Pe) +
    ylab(log10Po)
}
