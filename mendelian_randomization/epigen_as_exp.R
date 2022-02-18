library(TwoSampleMR)

library(tidyverse)

#Replace ieaa with hannum

exp <- read_exposure_data(
  "ieaa.txt",
  sep = "\t",
  beta_col = "Effect",
  snp_col = "rsID",
  se_col = "SE",
  effect_allele_col = "A1",
  other_allele_col = "A2",
  pval_col = "P",
  eaf_col = "Freq1",
  samplesize_col = "N"
)

#Replace ch with dnmt3a/tet2/large/small

out <- read_outcome_data(
  filename = "ch_mr_input.txt",
  snps = exp$SNP,
  sep = "\t",
  snp_col = "SNP",
  beta_col = "BETA1",
  se_col = "SE1",
  eaf_col = "A1FREQ",
  effect_allele_col = "ALLELE1",
  other_allele_col = "ALLELE0",
  pval_col = "P_BOLT_LMM",
)

dat <- harmonise_data(exp, out)

res <- generate_odds_ratios(mr(dat))

#Replace ch with dnmt3a/tet2/large/small
#Replace ieaa with hannum

write.table(res, file = "ieaa_as_exp_ch_out.txt", row.names = F, quote = F, sep = "\t")