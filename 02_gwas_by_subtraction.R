# Step 2: GWAS-by-subtraction of BMI from one OA trait.
# Usage: Rscript 02_gwas_by_subtraction.R <trait>    e.g. AllOA
#
# Output per trait:
#   RF    = OA genetic effects shared with BMI (BMI-related)
#   indep = OA genetic effects not mediated by BMI (BMI-subtracted)

library(GenomicSEM)
library(data.table)
library(dplyr)
library(readr)
source("config.R")

trait <- commandArgs(trailingOnly = TRUE)[1]
stopifnot(trait %in% traits)
out <- file.path(wd, "output", trait)
for (d in c("ldsc", "full_gwas", "rf_chunk", "indep_chunk", "RFlambda", "indeplambda"))
  dir.create(file.path(out, d), recursive = TRUE, showWarnings = FALSE)

# ---- LDSC on munged summary statistics ----
munged <- file.path(wd, "data/munge/munge_output", paste0(c(trait, "BMI"), ".sumstats.gz"))
setwd(file.path(out, "ldsc"))
LDSCoutput <- ldsc(munged, sample.prev = c(NA, NA), population.prev = c(NA, NA),
                   ld = ld, wld = ld, trait.names = c(trait, "BMI"))
save(LDSCoutput, file = file.path(out, "ldsc", paste0("LDSCoutput", trait, "-BMI.RData")))

# ---- Prepare full (not munged) summary statistics ----
to_sumstats_format <- function(df) {
  df %>%
    rename(SNP = rs_id, CHR = chromosome, POS = base_pair_location,
           A1 = effect_allele, A2 = other_allele, AF1 = effect_allele_frequency,
           BETA = beta, SE = standard_error, P = p_value, N = n) %>%
    select(SNP, CHR, POS, A1, A2, AF1, BETA, SE, P, N)
}
full <- file.path(out, "full_gwas", paste0(c(trait, "BMI"), "_cleaned.txt"))
write_tsv(to_sumstats_format(fread(oa_gwas(trait))), full[1])
write_tsv(to_sumstats_format(fread(bmi_gwas)), full[2])

p_sumstats <- sumstats(full, ref, trait.names = c(trait, "BMI"), maf.filter = 0.001,
                       se.logit = c(TRUE, FALSE), OLS = c(FALSE, TRUE), linprob = NULL,
                       N = c(LDSCoutput$N[1], LDSCoutput$N[3]))

# ---- GWAS-by-subtraction model ----
model <- gsub("OA", trait, "
  LatentOARF =~ NA*OA + BMI
  LatentOA   =~ NA*OA
  LatentOARF ~ SNP
  LatentOA   ~ SNP
  LatentOARF ~~ 1*LatentOARF
  LatentOA   ~~ 1*LatentOA
  LatentOA   ~~ 0*LatentOARF
  BMI ~~ 0*OA
  BMI ~~ 0*BMI
  OA  ~~ 0*OA
  SNP ~~ SNP")
sub <- gsub("OA", trait, c("LatentOARF~SNP", "LatentOA~SNP",
                           "LatentOARF=~OA", "LatentOARF=~BMI", "LatentOA=~OA"))

# Run in chunks; results are merged and rescaled in 03_merge_chunks.R
n_chunks <- ceiling(nrow(p_sumstats) / chunk_size)
for (chunk in seq_len(n_chunks)) {
  rows <- ((chunk - 1) * chunk_size + 1):min(chunk * chunk_size, nrow(p_sumstats))
  res <- userGWAS(covstruc = LDSCoutput, SNPs = p_sumstats[rows, ], estimation = "DWLS",
                  cores = cores, model = model, sub = sub)
  write_chunk <- function(x, d) fwrite(x, file.path(out, d, paste0("chunk", chunk, ".txt")), sep = "\t")
  write_chunk(res[[1]], "rf_chunk")     # LatentRF ~ SNP
  write_chunk(res[[2]], "indep_chunk")  # Latent ~ SNP
  write_chunk(res[[3]], "RFlambda")     # loading of OA on LatentRF
  write_chunk(res[[5]], "indeplambda")  # loading of OA on Latent
}
