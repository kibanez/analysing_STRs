# Objective: check how many genomes are being improved after QC visualisation
library(dplyr)

setwd("~/Documents/STRs/VALIDATION/QC_visual_inspection/raw_data_for_computing/")

l_68_ehv3 = read.table("./l_68_EHv3.txt", stringsAsFactors = F)
l_68_ehv3 = l_68_ehv3$V1
length(unique(l_68_ehv3))
# 65

# Fix `LP2000865-DNA_G07` since it contains ` `
# the same for # "LP3000543-DNA_C11" "LP3001453-DNA_B03" "LP3000449-DNA_C10"
all_table = read.csv("ehv3_before_after.tsv",
                     sep = "\t", 
                     stringsAsFactors = F, 
                     header = T)
dim(all_table)
# 634  14

only_checked = all_table %>% 
  filter(LP_number %in% l_68_ehv3)
dim(only_checked)
# 105  14

length(unique(only_checked$LP_number))
# 65

# How many genomes are before and after?
# TP
l_genomes_TP_before1 = only_checked %>% filter(before1 %in% "TP") %>% select(LP_number) %>% unique() %>% pull()
l_genomes_TP_before2 = only_checked %>% filter(before2 %in% "TP") %>% select(LP_number) %>% unique() %>% pull()
l_genomes_TP_before = unique(c(l_genomes_TP_before1, l_genomes_TP_before2))
length(l_genomes_TP_before)
# 55

# FP
l_genomes_FP_before1 = only_checked %>% filter(before1 %in% "FP") %>% select(LP_number) %>% unique() %>% pull()
l_genomes_FP_before2 = only_checked %>% filter(before2 %in% "FP") %>% select(LP_number) %>% unique() %>% pull()
l_genomes_FP_before = unique(c(l_genomes_FP_before1, l_genomes_FP_before2))
length(l_genomes_FP_before)
# 9

#TN
l_genomes_TN_before1 = only_checked %>% filter(before1 %in% "TN") %>% select(LP_number) %>% unique() %>% pull()
l_genomes_TN_before2 = only_checked %>% filter(before2 %in% "TN") %>% select(LP_number) %>% unique() %>% pull()
l_genomes_TN_before = unique(c(l_genomes_TN_before1, l_genomes_TN_before2))
length(l_genomes_TN_before)
# 54

# FN
l_genomes_FN_before1 = only_checked %>% filter(before1 %in% "FN") %>% select(LP_number) %>% unique() %>% pull()
l_genomes_FN_before2 = only_checked %>% filter(before2 %in% "FN") %>% select(LP_number) %>% unique() %>% pull()
l_genomes_FN_before = unique(c(l_genomes_FN_before1, l_genomes_FN_before2))
length(l_genomes_FN_before)
# 4

all_genomes = unique(c(l_genomes_TP_before,
                       l_genomes_FP_before,
                       l_genomes_FN_before))
length(all_genomes)
# 65

## Now after
# TP
l_genomes_TP_after1 = only_checked %>% filter(after1 %in% "TP") %>% select(LP_number) %>% unique() %>% pull()
l_genomes_TP_after2 = only_checked %>% filter(after2 %in% "TP") %>% select(LP_number) %>% unique() %>% pull()
l_genomes_TP_after = unique(c(l_genomes_TP_after1, l_genomes_TP_after2))
length(l_genomes_TP_after)
# 56

# FP
l_genomes_FP_after1 = only_checked %>% filter(after1 %in% "FP") %>% select(LP_number) %>% unique() %>% pull()
l_genomes_FP_after2 = only_checked %>% filter(after2 %in% "FP") %>% select(LP_number) %>% unique() %>% pull()
l_genomes_FP_after = unique(c(l_genomes_FP_after1, l_genomes_FP_after2))
length(l_genomes_FP_after)
# 2

#TN
l_genomes_TN_after1 = only_checked %>% filter(after1 %in% "TN") %>% select(LP_number) %>% unique() %>% pull()
l_genomes_TN_after2 = only_checked %>% filter(after2 %in% "TN") %>% select(LP_number) %>% unique() %>% pull()
l_genomes_TN_after = unique(c(l_genomes_TN_after1, l_genomes_TN_after2))
length(l_genomes_TN_after)
# 54

# FN
l_genomes_FN_after1 = only_checked %>% filter(after1 %in% "FN") %>% select(LP_number) %>% unique() %>% pull()
l_genomes_FN_after2 = only_checked %>% filter(after2 %in% "FN") %>% select(LP_number) %>% unique() %>% pull()
l_genomes_FN_after = unique(c(l_genomes_FN_after1, l_genomes_FN_after2))
length(l_genomes_FN_after)
# 0

all_genomes = unique(c(l_genomes_TP_after,
                       l_genomes_FP_after,
                       l_genomes_FN_after))
length(all_genomes)
# 57

# which 8 are missing?
setdiff(l_68_ehv3, all_genomes)
# "LP2000712-DNA_F10" "LP3000646-DNA_H06" "LP2000913-DNA_E03" "LP3000002-DNA_G02" "LP3000090-DNA_G06" "LP3000179-DNA_B02" "LP3000327-DNA_D01" "LP3000209-DNA_G02"
