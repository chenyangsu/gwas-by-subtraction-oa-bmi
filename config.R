# Edit these paths before running.

wd <- "/path/to/gwas-by-subtraction-oa-bmi"

# Input GWAS summary statistics in GWAS Catalog format, with columns:
# rs_id, chromosome, base_pair_location, effect_allele, other_allele,
# effect_allele_frequency, beta, standard_error, p_value, n
bmi_gwas <- "/path/to/BMI.tsv.gz"                    # Yengo et al. 2018
oa_gwas  <- function(trait) paste0("/path/to/oa_gwas/", trait, ".txt.gz")  # Boer et al. 2021

# Reference files from https://github.com/GenomicSEM/GenomicSEM/wiki
hm3 <- file.path(wd, "data/w_hm3.noMHC.snplist")
ld  <- file.path(wd, "data/eur_w_ld_chr/")
ref <- file.path(wd, "data/reference.1000G.maf.0.005.txt.gz")

# Traits and total sample sizes
traits <- c("AllOA", "KneeHipOA", "KneeOA", "HipOA", "TJR", "THR", "TKR",
            "SpineOA", "FingerOA", "ThumbOA", "HandOA", "Early_AllOA")
N <- c(BMI = 681275,
       AllOA = 826690, KneeHipOA = 490345, KneeOA = 396054, HipOA = 353388,
       TJR = 368576, THR = 319037, TKR = 252041, SpineOA = 333950,
       FingerOA = 266618, ThumbOA = 247455, HandOA = 303782, Early_AllOA = 48287)

chunk_size <- 10000
cores <- 10
