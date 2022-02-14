library(tidyverse)

library(vroom)

library(TwoSampleMR)

d <- vroom(file = "chTET2_0.6info_0.002MAF_filtered.txt.gz", delim = "\t")

d <- d %>% filter(INFO > 0.6)

d <- d %>% filter(A1FREQ > 0.01)

d <- d %>% filter(A1FREQ < 0.99)

#Change P_BOLT_LMM to < 5e-8 to use genome-wide significant variants only

d <- d %>% filter(P_BOLT_LMM < 1e-5)

d <- d %>% rename(chr_name = CHR, chrom_start = BP, pval.exposure = P_BOLT_LMM)

d <- clump_data(
  d,
  clump_kb = 10000,
  clump_r2 = 0.001,
  clump_p1 = 1,
  clump_p2 = 1,
  pop = "EUR"
)

d <- d %>% rename(BETA0 = BETA, SE0 = SE)

d <- d %>% mutate(BETA = BETA0/((2041/(2041+173918))*(1-(2041/(2041+173918)))))

d <- d %>% mutate(SE = SE0/((2041/(2041+173918))*(1-(2041/(2041+173918)))))

exp <- format_data(
  d,
  type = "exposure",
  snps = NULL,
  header = TRUE,
  snp_col = "SNP",
  beta_col = "BETA",
  se_col = "SE",
  eaf_col = "A1FREQ",
  effect_allele_col = "ALLELE1",
  other_allele_col = "ALLELE0",
  pval_col = "pval.exposure",
  chr_col = "chr_name",
  pos_col = "chrom_start"
)

#Download UKB_telomere_gwas_summarystats.tsv from https://figshare.com/s/caa99dc0f76d62990195

out <- read_outcome_data(
  filename = "UKB_telomere_gwas_summarystats.tsv",
  snps = exp$SNP,
  sep = "\t",
  snp_col = "variant_id",
  beta_col = "beta",
  se_col = "standard_error",
  eaf_col = "effect_allele_frequency",
  effect_allele_col = "effec_allele",
  other_allele_col = "other_allele",
  pval_col = "p_value"
)

dat <- harmonise_data(
  exposure_dat = exp,
  outcome_dat = out
)

res <- generate_odds_ratios(mr(dat))

write.table(res, file = "tet2_as_exp_tl_out.txt", row.names = F, quote = F, sep = "\t")