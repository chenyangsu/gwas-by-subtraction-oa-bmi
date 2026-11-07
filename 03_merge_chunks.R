# Step 3: merge chunks and rescale estimates to the scale of the OA trait.
# Usage: Rscript 03_merge_chunks.R <trait>    e.g. AllOA
#
# userGWAS returns SNP effects on standardized latent factors. Multiplying by the
# factor loading (lambda) puts them back on the OA scale:
#   est = est * lambda,  SE = SE * |lambda|

library(data.table)
source("config.R")

trait <- commandArgs(trailingOnly = TRUE)[1]
stopifnot(trait %in% traits)
out <- file.path(wd, "output", trait)

read_chunks <- function(d) {
  files <- list.files(file.path(out, d), pattern = "^chunk[0-9]+\\.txt$", full.names = TRUE)
  files <- files[order(as.integer(gsub("\\D", "", basename(files))))]
  rbindlist(lapply(files, fread))
}

rescale <- function(gwas, lambda) {
  gwas$est <- gwas$est * lambda$est
  gwas$SE  <- gwas$SE * abs(lambda$est)
  gwas$Z_Estimate <- gwas$est / gwas$SE
  gwas$Pval_Estimate <- 2 * pnorm(-abs(gwas$Z_Estimate))
  gwas
}

rf    <- rescale(read_chunks("rf_chunk"),    read_chunks("RFlambda"))
indep <- rescale(read_chunks("indep_chunk"), read_chunks("indeplambda"))

dir.create(file.path(out, "final_merged"), showWarnings = FALSE)
fwrite(rf,    file.path(out, "final_merged", paste0("outputGWAS_", trait, "_BMI_RF.txt")),    sep = "\t")
fwrite(indep, file.path(out, "final_merged", paste0("outputGWAS_", trait, "_BMI_indep.txt")), sep = "\t")
