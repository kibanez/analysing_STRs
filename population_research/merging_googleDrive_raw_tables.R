# Objective: from the recoded and formatted tables from google drive (https://docs.google.com/spreadsheets/d/1cuh2rsDkQP3YEHjX6ogLWlgxO3sd0Jfs4zCVdizBQEc/edit#gid=1383787997)
# The aim here is to dedup all info we do have from 3 independent tables and merge them all
date ()
Sys.info ()[c("nodename", "user")]
commandArgs ()
rm (list = ls ())
R.version.string ## "R version 3.6.1 (2019-07-05)"

library(ggplot2); packageDescription ("ggplot2", fields = "Version") #"2.2.1"
library(grid); packageDescription ("grid", fields = "Version") #"3.3.2"
library(gridExtra); packageDescription ("gridExtra", fields = "Version") #"2.3"
library(reshape2); packageDescription ("reshape2", fields = "Version") #"1.4.3"
require(dplyr); packageDescription ("dplyr", fields = "Version") #"0.7.4"

# Set working directory
setwd("~/Documents/STRs/VALIDATION/PCR_EH_estimations/")

# Load 3 tables
gel_table = read.csv("./googleDrive_GEL_validation.tsv",
                     sep = "\t",
                     stringsAsFactors = F, 
                     header = T)
dim(gel_table)
# 635 11

ari_table = read.csv("./googleDrive_Arianna_NHNN_validation.tsv",
                     sep = "\t",
                     stringsAsFactors = F, 
                     header = T)
dim(ari_table)
# 144 11 

james_table = read.csv("./googleDrive_James_NHNN_validation.tsv",
                      sep = "\t",
                      stringsAsFactors = F, 
                      header = T)
dim(james_table)
# 48  11

merge_all = rbind(gel_table,
                  ari_table,
                  james_table)
dim(merge_all)
# 827  11

merge_all = unique(merge_all)
dim(merge_all)
# 820  11

# which are duplicates?
intersect(gel_table$LP_number, ari_table$LP_number)
# "LP3000999-DNA_C06" "LP3001031-DNA_H09" "LP3000329-DNA_E12" "LP3001101-DNA_H06" "LP3000118-DNA_E03" "LP3000124-DNA_D08"
intersect(gel_table$LP_number, james_table$LP_number)
# "LP3000595-DNA_E05" "LP3000469-DNA_C05" "LP3000474-DNA_E01"


# remove all NA's, na's or '.'s
merge_all = merge_all %>%
  filter(!is.na(min.PCR.a1))
dim(merge_all)
# 810  11

merge_all = merge_all %>%
  filter(min.PCR.a1 != '.')
dim(merge_all)
# 792  11

# QC check -- they all have the same PCR< EHv2 and EHv3, estimations
write.table(merge_all, "googleDrive_all_merged_dedup_table.tsv", sep = "\t", quote = F, col.names = T, row.names = F)

# Analysis of the PCR-EH dataset
raw_numbers_popus = as.data.frame(table(merge_all$Super.population))
colnames(raw_numbers_popus) = c("population", "Number of genomes")


png("figures/barplot_ancestry_PCR_cohort.png")
ggplot(raw_numbers_popus, 
       aes(x = reorder(population, -`Number of genomes`), y = `Number of genomes`)) + 
  geom_bar(stat = "identity", aes(fill = population)) + 
  geom_text(aes(label=`Number of genomes`), vjust=-0.5, size = 4, colour = "grey") +
  ylab("Number of genomes") + 
  xlab("Available PCR cohort - 792 genomes") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
dev.off()

# Number of genomes per locus
raw_numbers_locus = as.data.frame(table(merge_all$locus))
colnames(raw_numbers_locus) = c("locus", "Number of genomes")

png("figures/barplot_loci_distribution_PCR_cohort.png")
ggplot(raw_numbers_locus, 
       aes(x = reorder(locus, -`Number of genomes`), y = `Number of genomes`)) + 
  geom_bar(stat = "identity", aes(fill = locus)) + 
  geom_text(aes(label=`Number of genomes`), vjust=-0.5, size = 4, colour = "grey") +
  ylab("Number of genomes") + 
  xlab("Available PCR cohort - 792 genomes") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
dev.off()

# Bubble plot with PCR estimations vs EHv255

# Prepare data per each superpopu
l_superpopu = unique(merge_all$Super.population)
l_superpopu = l_superpopu[-5]
for (j in 1:length(l_superpopu)){
 superpopu_merge = merge_all %>% filter(Super.population %in% l_superpopu[j])
  
  # Create dataframe with exp, eh, freq for each locus
  df_data_with_freq_v2 = data.frame()
  l_locus = unique(superpopu_merge$locus)
  for(i in 1:length(l_locus)){
    aux_validation_a1 = superpopu_merge %>% filter(locus %in% l_locus[i]) %>% select(min.PCR.a1) %>% pull() %>% as.integer() 
    aux_validation_a2 = superpopu_merge %>% filter(locus %in% l_locus[i]) %>% select(maxPCR.a2) %>% pull() %>% as.integer() 
    aux_exp_alleles_v2 = c(aux_validation_a1, aux_validation_a2)
    
    aux_eh_a1 = superpopu_merge %>% filter(locus %in% l_locus[i]) %>% select(min.EHv255.a1) %>% pull() %>% as.integer() 
    aux_eh_a2 = superpopu_merge %>% filter(locus %in% l_locus[i]) %>% select(max.EHv255.a2) %>% pull() %>% as.integer() 
    aux_eh_alleles_v2 = c(aux_eh_a1, aux_eh_a2)
    
    data_aux = xyTable(aux_exp_alleles_v2, aux_eh_alleles_v2)
    
    df_data_aux = data.frame(eh_alleles = data_aux$y,
                             exp_alleles = data_aux$x,
                             number_of_alleles = data_aux$number,
                             locus = rep(l_locus[i], length(data_aux$x)))
    # Concat info per locus
    df_data_with_freq_v2 = rbind(df_data_with_freq_v2,
                                 df_data_aux)
    
  }
  
  max_value = max(df_data_with_freq_v2$eh_alleles, 
                  df_data_with_freq_v2$exp_alleles,
                  na.rm = TRUE) + 5
  
  group.colors = c("AR" = "#1B9E77", "ATN1" = "#D95F02", "ATXN1" ="#7570B3", "ATXN2" = "#E7298A", "ATXN3" = "#66A61E", 
                   "ATXN7" = "#E6AB02", "CACNA1A" = "#A6761D", "FXN" = "#666666", "HTT" ="#E41A1C", "TBP" = "#FF7F00", 
                   "C9orf72" = "#FFFF33", "FMR1" = "#F781BF", "PPP2R2B" = "black")
  
  
  joint_plot = ggplot(df_data_with_freq_v2, 
                      aes(x = exp_alleles, y = eh_alleles, colour = factor(locus))) + 
    geom_point(aes(fill = factor(locus), size = number_of_alleles)) + 
    xlim(5,max_value) + 
    ylim(5,max_value) + 
    labs(title = "", 
         y = "EHv255 repeat sizes", 
         x = "PCR repeat sizes") + 
    geom_abline(method = "lm", formula = x ~ y, linetype = 2, colour = "gray") +  
    coord_equal() +
    scale_fill_manual(values=group.colors) +  
    theme(legend.title = element_blank(),
          axis.text.x.top = element_text()) + 
    guides(size = FALSE) 
  
  
  figure_name = paste("./figures/Figure1_", l_superpopu[j], sep = "")
  figure_name = paste(figure_name, "across_loci.png", sep = "_")
  png(figure_name,units="in", width=5, height=5, res=600)
  print(joint_plot)
  dev.off()
  
  
  # breakdown by locus
  joint_plot_breakdown = ggplot(df_data_with_freq_v2, 
                      aes(x = exp_alleles, y = eh_alleles, colour = factor(locus))) + 
    geom_point(aes(fill = factor(locus), size = number_of_alleles)) + 
    xlim(5,max_value) + 
    ylim(5,max_value) + 
    labs(title = "", 
         y = "EHv255 repeat sizes", 
         x = "PCR repeat sizes") + 
    geom_abline(method = "lm", formula = x ~ y, linetype = 2, colour = "gray") +  
    coord_equal() +
    scale_fill_manual(values=group.colors) +  
    theme(legend.title = element_blank(),
          axis.text.x.top = element_text()) + 
    guides(size = FALSE) +
    facet_wrap(locus~ .)

  figure_name = paste("./figures/Figure1_", l_superpopu[j], sep = "")
  figure_name = paste(figure_name, "across_loci_brokendown_by_loci.png", sep = "_")
  png(figure_name,units="in", width=5, height=5, res=600)
  print(joint_plot_breakdown)
  dev.off()
  
  
}# superpopu

# Bubble plot with PCR estimations vs EHv312