# GWAS-by-subtraction: osteoarthritis minus BMI

Code for the GWAS-by-subtraction analysis in:

> Su C-Y, et al. *Disentangling osteoarthritis-specific genetic effects from obesity to identify novel therapeutic targets.* Nature Communications.

Each osteoarthritis (OA) GWAS is split into two components with
[GenomicSEM](https://github.com/GenomicSEM/GenomicSEM) (Demange et al., 2021):

- **RF (BMI-related):** OA genetic effects shared with BMI
- **indep (BMI-subtracted):** OA genetic effects not mediated by BMI

## Usage

1. Edit the paths in `config.R`.
2. Run:

```bash
Rscript 01_munge.R
for t in AllOA KneeHipOA KneeOA HipOA TJR THR TKR SpineOA FingerOA ThumbOA HandOA Early_AllOA; do
  Rscript 02_gwas_by_subtraction.R $t
  Rscript 03_merge_chunks.R $t
done
```

Results are written to `output/<trait>/final_merged/outputGWAS_<trait>_BMI_{RF,indep}.txt`.

| Script | What it does |
|---|---|
| `config.R` | Paths, traits, sample sizes |
| `01_munge.R` | Munges BMI and OA summary statistics to HapMap3 SNPs (INFO > 0.9, MAF > 0.001) |
| `02_gwas_by_subtraction.R` | LDSC, `sumstats`, and the GWAS-by-subtraction model with `userGWAS`, in chunks |
| `03_merge_chunks.R` | Merges chunks and rescales estimates and SEs by the factor loadings |

## Model

```
LatentOARF =~ NA*OA + BMI
LatentOA   =~ NA*OA
LatentOARF ~ SNP
LatentOA   ~ SNP
LatentOARF ~~ 1*LatentOARF
LatentOA   ~~ 1*LatentOA
LatentOA   ~~ 0*LatentOARF
BMI ~~ 0*OA;  BMI ~~ 0*BMI;  OA ~~ 0*OA;  SNP ~~ SNP
```

## Requirements

R (tested with 4.1.2) with GenomicSEM, data.table, dplyr and readr.

Reference files (`w_hm3.noMHC.snplist`, `eur_w_ld_chr/`, `reference.1000G.maf.0.005.txt.gz`)
are available from the [GenomicSEM wiki](https://github.com/GenomicSEM/GenomicSEM/wiki).

Input GWAS: OA from the Genetics of Osteoarthritis Consortium (Boer et al., 2021) and
BMI from GIANT-UK Biobank (Yengo et al., 2018).

## Data

The 24 resulting summary statistics are deposited in the NHGRI-EBI GWAS Catalog
(GCST91006158 to GCST91006181).
