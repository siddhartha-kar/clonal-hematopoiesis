library(TwoSampleMR)

library(tidyverse)

library(vroom)

ids <- vroom(file = "ahola-olli_ids.txt", delim = "\t")

exp <- extract_instruments(outcomes=ids$id, p1 = 1e-05, clump = TRUE, p2 = 1e-05, r2 = 0.001, kb = 10000)

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

res <- res %>% arrange(method, pval)

#Replace ch with dnmt3a/tet2/large/small

write.table(res, file = "cytokines_gfs_as_exp_ch_out.txt", row.names = F, quote = F, sep = "\t")