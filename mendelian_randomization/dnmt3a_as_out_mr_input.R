library(tidyverse)

library(vroom)

d <- vroom(file = "chipDNMT3A_bolt_0.3info_0.00005MAF_filtered.out", delim = "\t")

d <- d %>% filter(INFO > 0.3)

d <- d %>% filter(A1FREQ > 0.01)

d <- d %>% filter(A1FREQ < 0.99)

d <- d %>% mutate(BETA1 = BETA/((5185/(5185+173918))*(1-(5185/(5185+173918)))))

d <- d %>% mutate(SE1 = SE/((5185/(5185+173918))*(1-(5185/(5185+173918)))))

write.table(d, file = "dnmt3a_mr_input.txt", sep = "\t", row.names = F, quote = F)