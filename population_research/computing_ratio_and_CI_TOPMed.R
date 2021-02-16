# Objective: from TOPMed tables, compute the carrier freq ratio and CI values
date ()
Sys.info ()[c("nodename", "user")]
commandArgs ()
rm (list = ls ())
R.version.string ## "R version 3.6.1 (2019-07-05)"

# libraries
library(dplyr)

# Set working directory
setwd("~/Documents/STRs/ANALYSIS/population_research/TOPMed/After_QC/")

total_unrel = 49239
total_unrel_notNeuro = 48389

# Load table
topmed_unrel = read.csv("./table_unrel_and_unrel_notNeuro.tsv", stringsAsFactors = F, header = T, sep = "\t")
dim(topmed_unrel)
# 13  3

l_locus = unique(topmed_unrel$locus)
df_unrel = data.frame()
df_unrel_notNeuro = data.frame()
for(i in 1:length(l_locus)){
  # Unrel
  # Compute number of expanded genomes per locus (after visual inspection)
  total_exp_after_VI_locus = topmed_unrel %>%
    filter(locus %in% l_locus[i]) %>%
    select(exp_after_VI_unrel) %>%
    unique() %>%
    pull() 
  
  freq_carrier = round(total_unrel / total_exp_after_VI_locus,digits = 2)
  ratio_freq_carrier = paste("1 in", as.character(freq_carrier), sep = " ")
  
  ci_max = round(total_unrel/(total_unrel*((total_exp_after_VI_locus/total_unrel)-1.96*sqrt((total_exp_after_VI_locus/total_unrel)*(1-total_exp_after_VI_locus/total_unrel)/total_unrel))), digits = 2)
  ci_min = round(total_unrel/(total_unrel*((total_exp_after_VI_locus/total_unrel)+1.96*sqrt((total_exp_after_VI_locus/total_unrel)*(1-total_exp_after_VI_locus/total_unrel)/total_unrel))), digits = 2)
  
  ci_ratio= as.character(paste(as.character(ci_min), as.character(ci_max), sep = "-"))
  
  df_unrel = rbind(df_unrel,
                   cbind(l_locus[i], total_exp_after_VI_locus, total_unrel, ratio_freq_carrier, ci_ratio))
  
  # Unrel NOT NEURO
  # Compute number of expanded genomes per locus (after visual inspection)
  total_exp_after_VI_locus_notNeuro = topmed_unrel %>%
    filter(locus %in% l_locus[i]) %>%
    select(exp_after_VI_unrel_notNeuro) %>%
    unique() %>%
    pull()
  
  freq_carrier_notNeuro = round(total_unrel_notNeuro / total_exp_after_VI_locus_notNeuro,digits = 2)
  ratio_freq_carrier_notNeuro = paste("1 in", as.character(freq_carrier_notNeuro), sep = " ")
  
  ci_max_notNeuro = round(total_unrel_notNeuro/(total_unrel_notNeuro*((total_exp_after_VI_locus_notNeuro/total_unrel_notNeuro)-1.96*sqrt((total_exp_after_VI_locus_notNeuro/total_unrel_notNeuro)*(1-total_exp_after_VI_locus_notNeuro/total_unrel_notNeuro)/total_unrel_notNeuro))), digits = 2)
  ci_min_notNeuro = round(total_unrel_notNeuro/(total_unrel_notNeuro*((total_exp_after_VI_locus_notNeuro/total_unrel_notNeuro)+1.96*sqrt((total_exp_after_VI_locus_notNeuro/total_unrel_notNeuro)*(1-total_exp_after_VI_locus_notNeuro/total_unrel_notNeuro)/total_unrel_notNeuro))), digits = 2)
  
  ci_ratio_notNeuro = as.character(paste(as.character(ci_min_notNeuro), as.character(ci_max_notNeuro), sep = "-"))
  
  df_unrel_notNeuro = rbind(df_unrel_notNeuro,
                            cbind(l_locus[i], total_exp_after_VI_locus_notNeuro, total_unrel_notNeuro, ratio_freq_carrier_notNeuro, ci_ratio_notNeuro))
  
}

# write them into files
write.table(df_unrel,
            "table_carrier_ratio_with_ci_unrel_genomes_TOPMed.tsv",
            quote = F,
            col.names = T,
            row.names = F,
            sep = "\t")

write.table(df_unrel_notNeuro,
            "table_carrier_ratio_with_ci_unrel_NotNeuro_genomes_TOPMed.tsv",
            quote = F,
            col.names = T,
            row.names = F,
            sep = "\t")








