library(TwoSampleMR)

library(tidyverse)

library(vroom)

ids <- vroom(file = "astle_ids.txt", delim = "\t")

exp <- extract_instruments(outcomes=ids$id)

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

write.table(res, file = "blood_traits_as_exp_ch_out.txt", row.names = F, quote = F, sep = "\t")