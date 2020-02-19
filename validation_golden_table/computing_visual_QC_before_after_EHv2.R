# Objective: check how many genomes are being improved after QC visualisation
library(dplyr)

setwd("~/Documents/STRs/VALIDATION/QC_visual_inspection/raw_data_for_computing/")

l_66_ehv2 = read.table("list_66_EHv2.txt", stringsAsFactors = F)
l_66_ehv2 = l_66_ehv2$V1
length(unique(l_66_ehv2))
# 66

# Fix `LP2000865-DNA_G07` since it contains ` `
all_table = read.csv("~/Downloads/ehv2_before_after.tsv",
                     sep = "\t", 
                     stringsAsFactors = F, 
                     header = T)
dim(all_table)
# 634  14

only_checked = all_table %>% 
  filter(LP_number %in% l_66_ehv2)
dim(only_checked)
# 120  14

length(unique(only_checked$LP_number))
# 66

# Which 3 genomes are missing here??
setdiff(l_66_ehv2,unique(only_checked$LP_number))
# "LP3000543-DNA_C11" "LP3001453-DNA_B03" "LP3000449-DNA_C10"


# How many genomes are before and after?
# TP
l_genomes_TP_before1 = only_checked %>% filter(before1 %in% "TP") %>% select(LP_number) %>% unique() %>% pull()
l_genomes_TP_before2 = only_checked %>% filter(before2 %in% "TP") %>% select(LP_number) %>% unique() %>% pull()
l_genomes_TP_before = unique(c(l_genomes_TP_before1, l_genomes_TP_before2))
length(l_genomes_TP_before)
# 56

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
# 55

# FN
l_genomes_FN_before1 = only_checked %>% filter(before1 %in% "FN") %>% select(LP_number) %>% unique() %>% pull()
l_genomes_FN_before2 = only_checked %>% filter(before2 %in% "FN") %>% select(LP_number) %>% unique() %>% pull()
l_genomes_FN_before = unique(c(l_genomes_FN_before1, l_genomes_FN_before2))
length(l_genomes_FN_before)
# 1

all_genomes = unique(c(l_genomes_TP_before,
                       l_genomes_FP_before,
                       l_genomes_FN_before))
length(all_genomes)
# 66

## Now after