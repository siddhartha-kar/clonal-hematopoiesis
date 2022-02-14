library(tidyverse)

library(vroom)

library(TwoSampleMR)

d <- vroom(file = "chOVERALL_0.6info_0.002MAF_filtered.txt.gz", delim = "\t")

d <- d %>% filter(INFO > 0.6)

d <- d %>% filter(A1FREQ > 0.01)

d <- d %>% filter(A1FREQ < 0.99)

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

d <- d %>% mutate(BETA = BETA0/((10203/(10203+173918))*(1-(10203/(10203+173918)))))

d <- d %>% mutate(SE = SE0/((10203/(10203+173918))*(1-(10203/(10203+173918)))))

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

ids <- vroom(file = "ukb-ad_outcomes_ids.txt", delim = "\t")

out <- extract_outcome_data(
  snps = exp$SNP,
  outcomes = ids$id
)

dat <- harmonise_data(
  exposure_dat = exp,
  outcome_dat = out
)

res <- generate_odds_ratios(mr(dat))

write.table(res, file = "ch_as_exp_ukb-ad_out.txt", row.names = F, quote = F, sep = "\t")
