# Step 1: munge BMI and OA summary statistics to HapMap3 SNPs for LDSC.
# Usage: Rscript 01_munge.R

library(GenomicSEM)
library(data.table)
library(dplyr)
library(readr)
source("config.R")

input_dir  <- file.path(wd, "data/munge/munge_input")
output_dir <- file.path(wd, "data/munge/munge_output")
dir.create(input_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

to_munge_format <- function(df) {
  df %>%
    rename(MarkerName = rs_id, CHR = chromosome, POS = base_pair_location,
           A1 = effect_allele, A2 = other_allele, EAF = effect_allele_frequency,
           Beta = beta, SE = standard_error, Pval = p_value) %>%
    select(MarkerName, CHR, POS, A1, A2, EAF, Beta, SE, Pval)
}

all_traits <- c("BMI", traits)
gwas_files <- c(bmi_gwas, sapply(traits, oa_gwas))
files <- file.path(input_dir, paste0(all_traits, ".txt"))
for (i in seq_along(all_traits)) {
  write_tsv(to_munge_format(fread(gwas_files[i])), files[i])
}

setwd(output_dir)  # munge() writes <trait>.sumstats.gz to the working directory
munge(files = files, hm3 = hm3, trait.names = all_traits, N = N[all_traits],
      info.filter = 0.9, maf.filter = 0.001)
