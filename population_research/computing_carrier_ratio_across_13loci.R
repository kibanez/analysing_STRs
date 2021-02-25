# Objective: from the work we have done by inspecting visually all pileups, compute the carrier ratio for each locus
# unrelated, unrelated not neuro
date ()
Sys.info ()[c("nodename", "user")]
commandArgs ()
rm (list = ls ())
R.version.string ## "R version 3.6.1 (2019-07-05)"

# libraries
library(dplyr)
library(tidyverse)

# Set working directory
setwd("~/Documents/STRs/ANALYSIS/population_research/PAPER/carriers/cc_pileup_100Kg/")

# Load unrel 55603 clin data
clin_data = read.csv("../table_55603_unrel_genomes_enriched_popu_diseasegroup.tsv",
                     stringsAsFactors = F,
                     header = F,
                     sep = "\t")
dim(clin_data)
# 55603  5

colnames(clin_data) = c("platekey", "famID", "disease_group", "is_neuro", "popu")

# summarise here total number of unrel genomes
# summarise here total number of unrel genomes for each enthnicity
total_unrel = length(unique(clin_data$platekey))
# 55603

# Load 2x125bp sequenced genomes
l_125 = read.table("../../list_genomes_125bp_100kGP.tsv", stringsAsFactors = F)
l_125 = l_125$V1
length(l_125)
# 15830

total_unrel = setdiff(unique(clin_data$platekey), l_125)
length(total_unrel)
# 54437 

total_unrel_AFR = clin_data %>% filter(popu %in% "AFR") %>% select(platekey) %>% unique() %>% pull() %>% length()
total_unrel_AMR = clin_data %>% filter(popu %in% "AMR") %>% select(platekey) %>% unique() %>% pull() %>% length()
total_unrel_EAS = clin_data %>% filter(popu %in% "EAS") %>% select(platekey) %>% unique() %>% pull() %>% length()
total_unrel_EUR = clin_data %>% filter(popu %in% "EUR") %>% select(platekey) %>% unique() %>% pull() %>% length()
total_unrel_SAS = clin_data %>% filter(popu %in% "SAS") %>% select(platekey) %>% unique() %>% pull() %>% length()

total_unrel_notNeuro = clin_data %>% filter(is_neuro %in% "NotNeuro") %>% select(platekey) %>% unique() %>% pull() %>% length()
# 37888
total_unrel_AFR_notNeuro = clin_data %>% filter(popu %in% "AFR", is_neuro %in% "NotNeuro") %>% select(platekey) %>% unique() %>% pull() %>% length()
total_unrel_AMR_notNeuro = clin_data %>% filter(popu %in% "AMR", is_neuro %in% "NotNeuro") %>% select(platekey) %>% unique() %>% pull() %>% length()
total_unrel_EAS_notNeuro = clin_data %>% filter(popu %in% "EAS", is_neuro %in% "NotNeuro") %>% select(platekey) %>% unique() %>% pull() %>% length()
total_unrel_EUR_notNeuro = clin_data %>% filter(popu %in% "EUR", is_neuro %in% "NotNeuro") %>% select(platekey) %>% unique() %>% pull() %>% length()
total_unrel_SAS_notNeuro = clin_data %>% filter(popu %in% "SAS", is_neuro %in% "NotNeuro") %>% select(platekey) %>% unique() %>% pull() %>% length()

# Load the whole table for 100kGP - case-controls 
table_100cc_QC = read.csv("./table_platekey_locus_QC_inspection_19feb21.tsv",
                          stringsAsFactors = F,
                          header = T,
                          sep = "\t")
dim(table_100cc_QC)
# 1783  16

# For each locus, compute the carrier ratio and CI 
l_locus = unique(table_100cc_QC$locus)
df_unrel = data.frame()
df_unrel_notNeuro = data.frame()
for(i in 1:length(l_locus)){
  # Unrel
  # Compute number of expanded genomes per locus (after visual inspection)
  total_exp_after_VI_locus = table_100cc_QC %>%
    filter(locus %in% l_locus[i], Final.decision %in% "Yes", is_unrelated. %in% "Yes") %>%
    select(platekey) %>%
    unique() %>%
    pull() %>%
    length()
  
  freq_carrier = round(total_unrel / total_exp_after_VI_locus,digits = 2)
  ratio_freq_carrier = paste("1 in", as.character(freq_carrier), sep = " ")
  
  ci_max = round(total_unrel/(total_unrel*((total_exp_after_VI_locus/total_unrel)-1.96*sqrt((total_exp_after_VI_locus/total_unrel)*(1-total_exp_after_VI_locus/total_unrel)/total_unrel))), digits = 2)
  ci_min = round(total_unrel/(total_unrel*((total_exp_after_VI_locus/total_unrel)+1.96*sqrt((total_exp_after_VI_locus/total_unrel)*(1-total_exp_after_VI_locus/total_unrel)/total_unrel))), digits = 2)
  
  ci_ratio= as.character(paste(as.character(ci_min), as.character(ci_max), sep = "-"))
  
  df_unrel = rbind(df_unrel,
                   cbind(l_locus[i], total_exp_after_VI_locus, total_unrel, ratio_freq_carrier, ci_ratio))

  # Unrel NOT NEURO
  # Compute number of expanded genomes per locus (after visual inspection)
  total_exp_after_VI_locus_notNeuro = table_100cc_QC %>%
    filter(locus %in% l_locus[i], Final.decision %in% "Yes", is_unrelated. %in% "Yes", is_neuro. %in% "NotNeuro") %>%
    select(platekey) %>%
    unique() %>%
    pull() %>%
    length()
  
  freq_carrier_notNeuro = round(total_unrel_notNeuro / total_exp_after_VI_locus_notNeuro,digits = 2)
  ratio_freq_carrier_notNeuro = paste("1 in", as.character(freq_carrier_notNeuro), sep = " ")
  
  ci_max_notNeuro = round(total_unrel_notNeuro/(total_unrel_notNeuro*((total_exp_after_VI_locus_notNeuro/total_unrel_notNeuro)-1.96*sqrt((total_exp_after_VI_locus_notNeuro/total_unrel_notNeuro)*(1-total_exp_after_VI_locus_notNeuro/total_unrel_notNeuro)/total_unrel_notNeuro))), digits = 2)
  ci_min_notNeuro = round(total_unrel_notNeuro/(total_unrel_notNeuro*((total_exp_after_VI_locus_notNeuro/total_unrel_notNeuro)+1.96*sqrt((total_exp_after_VI_locus_notNeuro/total_unrel_notNeuro)*(1-total_exp_after_VI_locus_notNeuro/total_unrel_notNeuro)/total_unrel_notNeuro))), digits = 2)
  
  ci_ratio_notNeuro = as.character(paste(as.character(ci_min_notNeuro), as.character(ci_max_notNeuro), sep = "-"))
  
  df_unrel_notNeuro = rbind(df_unrel_notNeuro,
                            cbind(l_locus[i], total_exp_after_VI_locus_notNeuro, total_unrel_notNeuro, ratio_freq_carrier_notNeuro, ci_ratio_notNeuro))
  
}

# Write unrel and unrel not neuro into files
write.table(df_unrel,
            "table_carrier_ratio_with_ci_unrel_genomes.tsv",
            quote = F,
            col.names = T,
            row.names = F,
            sep = "\t")

write.table(df_unrel_notNeuro,
            "table_carrier_ratio_with_ci_unrel_NotNeuro_genomes.tsv",
            quote = F,
            col.names = T,
            row.names = F,
            sep = "\t")

# By population
l_popu = c("AFR", "AMR", "EAS", "EUR", "SAS")
total_unrel_popu = c(total_unrel_AFR, total_unrel_AMR, total_unrel_EAS, total_unrel_EUR, total_unrel_SAS)
total_unrel_popu_notNeuro = c(total_unrel_AFR_notNeuro, total_unrel_AMR_notNeuro, total_unrel_EAS_notNeuro, total_unrel_EUR_notNeuro, total_unrel_SAS_notNeuro)
df_unrel = data.frame()
df_unrel_notNeuro = data.frame()
for(j in 1:length(l_popu)){
  df_unrel = data.frame()
  df_unrel_notNeuro = data.frame()
  for(i in 1:length(l_locus)){
  # Unrel
  # Compute number of expanded genomes per locus (after visual inspection)
  total_exp_after_VI_locus = table_100cc_QC %>%
    filter(locus %in% l_locus[i], Final.decision %in% "Yes", is_unrelated. %in% "Yes", popu %in% l_popu[j]) %>%
    select(platekey) %>%
    unique() %>%
    pull() %>%
    length()
  
  freq_carrier = round(total_unrel_popu[j] / total_exp_after_VI_locus,digits = 2)
  ratio_freq_carrier = paste("1 in", as.character(freq_carrier), sep = " ")
  
  ci_max = round(total_unrel_popu[j]/(total_unrel_popu[j]*((total_exp_after_VI_locus/total_unrel_popu[j])-1.96*sqrt((total_exp_after_VI_locus/total_unrel_popu[j])*(1-total_exp_after_VI_locus/total_unrel_popu[j])/total_unrel_popu[j]))), digits = 2)
  ci_min = round(total_unrel_popu[j]/(total_unrel_popu[j]*((total_exp_after_VI_locus/total_unrel_popu[j])+1.96*sqrt((total_exp_after_VI_locus/total_unrel_popu[j])*(1-total_exp_after_VI_locus/total_unrel_popu[j])/total_unrel_popu[j]))), digits = 2)
  
  ci_ratio= as.character(paste(as.character(ci_min), as.character(ci_max), sep = "-"))
  
  df_unrel = rbind(df_unrel,
                   cbind(l_locus[i], total_exp_after_VI_locus, total_unrel_popu[j], ratio_freq_carrier, ci_ratio))
  
  # Unrel NOT NEURO
  # Compute number of expanded genomes per locus (after visual inspection)
  total_exp_after_VI_locus_notNeuro = table_100cc_QC %>%
    filter(locus %in% l_locus[i], Final.decision %in% "Yes", is_unrelated. %in% "Yes", is_neuro. %in% "NotNeuro", popu %in% l_popu[j]) %>%
    select(platekey) %>%
    unique() %>%
    pull() %>%
    length()
  
  freq_carrier_notNeuro = round(total_unrel_popu_notNeuro[j] / total_exp_after_VI_locus_notNeuro,digits = 2)
  ratio_freq_carrier_notNeuro = paste("1 in", as.character(freq_carrier_notNeuro), sep = " ")
  
  ci_max_notNeuro = round(total_unrel_popu_notNeuro[j]/(total_unrel_popu_notNeuro[j]*((total_exp_after_VI_locus_notNeuro/total_unrel_popu_notNeuro[j])-1.96*sqrt((total_exp_after_VI_locus_notNeuro/total_unrel_popu_notNeuro[j])*(1-total_exp_after_VI_locus_notNeuro/total_unrel_popu_notNeuro[j])/total_unrel_popu_notNeuro[j]))), digits = 2)
  ci_min_notNeuro = round(total_unrel_popu_notNeuro[j]/(total_unrel_popu_notNeuro[j]*((total_exp_after_VI_locus_notNeuro/total_unrel_popu_notNeuro[j])+1.96*sqrt((total_exp_after_VI_locus_notNeuro/total_unrel_popu_notNeuro[j])*(1-total_exp_after_VI_locus_notNeuro/total_unrel_popu_notNeuro[j])/total_unrel_popu_notNeuro[j]))), digits = 2)
  
  ci_ratio_notNeuro = as.character(paste(as.character(ci_min_notNeuro), as.character(ci_max_notNeuro), sep = "-"))
  
  df_unrel_notNeuro = rbind(df_unrel_notNeuro,
                            cbind(l_locus[i], total_exp_after_VI_locus_notNeuro, total_unrel_popu_notNeuro[j], ratio_freq_carrier_notNeuro, ci_ratio_notNeuro))

  }
  # Write into files  - per popu
  output_file = paste(paste("table_carrier_ratio_with_ci_unrel_genomes", l_popu[j], sep = "_"), ".tsv", sep = "")
  write.table(df_unrel, output_file, quote = F, row.names = F, col.names = T, sep = "\t")
  output_file = paste(paste("table_carrier_ratio_with_ci_unrel__NotNeuro_genomes", l_popu[j], sep = "_"), ".tsv", sep = "")
  write.table(df_unrel_notNeuro, output_file, quote = F, row.names = F, col.names = T, sep = "\t")
}
