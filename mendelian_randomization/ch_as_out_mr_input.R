library(tidyverse)

library(vroom)

d <- vroom(file = "chip_bolt_0.3info_0.00005MAF_filtered.out", delim = "\t")

d <- d %>% filter(INFO > 0.3)

d <- d %>% filter(A1FREQ > 0.01)

d <- d %>% filter(A1FREQ < 0.99)

d <- d %>% mutate(BETA1 = BETA/((10203/(10203+173918))*(1-(10203/(10203+173918)))))

d <- d %>% mutate(SE1 = SE/((10203/(10203+173918))*(1-(10203/(10203+173918)))))

write.table(d, file = "ch_mr_input.txt", sep = "\t", row.names = F, quote = F)