# Objective: represent with bubble plots the correspondance between PCR and EH sizes
# split by super-population
date ()
Sys.info ()[c("nodename", "user")]
commandArgs ()
rm (list = ls ())
R.version.string ## "R version 3.6.3(2020-02-29)"

library(ggplot2); packageDescription ("ggplot2", fields = "Version") #"2.2.1"
library(grid); packageDescription ("grid", fields = "Version") #"3.3.2"
library(gridExtra); packageDescription ("gridExtra", fields = "Version") #"2.3"
library(reshape2); packageDescription ("reshape2", fields = "Version") #"1.4.3"
require(dplyr); packageDescription ("dplyr", fields = "Version") #"0.7.4"
library(lattice); packageDescription ("lattice", fields = "Version") #"0.20-35"
library(tidyverse); packageDescription ("tidyverse", fields = "Version") #"1.2.1"
library(RColorBrewer); packageDescription ("RColorBrewer", fields = "Version") #"1.1-2"
library(cowplot); packageDescription ("cowplot", fields = "Version") #"1.0.0"

# Set working environment
setwd("/Users/kibanez/Documents/STRs/PAPERS/POPULATION/figures/")

# Loading the golden or truthset or performance dataset
# Removing from here large expansions
# Only those < read-length
val_data = read.csv("~/Documents/STRs/ANALYSIS/population_research/PAPER/PerformanceDataset/PCR_vs_EH_all_together_POPU_PAPER.tsv", sep = "\t", stringsAsFactors = F, header = T)
dim(val_data)
# 714  17

l_superpopus = unique(val_data$superpopu)
for(i in 1:length(l_superpopus)){
  val_data_superpopu = val_data %>%
    filter(superpopu %in% l_superpopus[i])
  
  # Let's take the important meat
  exp_alleles_v2 = c(as.integer(val_data_superpopu$min.PCR.a1), as.integer(val_data_superpopu$maxPCR.a2))
  eh_alleles_v2 = c(as.integer(val_data_superpopu$min.EHv312.a1), as.integer(val_data_superpopu$max.EHv312.a2))
  locus_v2 = c(val_data_superpopu$locus, val_data_superpopu$locus)
  
  # Remove NAs
  index_NA = which(is.na(exp_alleles_v2))
  # Remove EXP,NORM,PREMUT,FULL EXP
  index_EXP = which(grepl("EXP",exp_alleles_v2))
  
  if (length(index_NA) > 0){
    exp_alleles_v2 = exp_alleles_v2[-index_NA]
    eh_alleles_v2 = eh_alleles_v2[-index_NA]
    locus_v2 = locus_v2[-index_NA]
    
  }
  
  # Create dataframe with exp, eh, freq for each locus
  df_data_with_freq_v2 = data.frame()
  l_locus = unique(locus_v2)
  for(j in 1:length(l_locus)){
    aux_validation_a1 = val_data_superpopu %>% filter(locus %in% l_locus[j]) %>% select(min.PCR.a1) %>% pull() %>% as.integer() 
    aux_validation_a2 = val_data_superpopu %>% filter(locus %in% l_locus[j]) %>% select(maxPCR.a2) %>% pull() %>% as.integer() 
    aux_exp_alleles_v2 = c(aux_validation_a1, aux_validation_a2)
    
    aux_eh_a1 = val_data_superpopu %>% filter(locus %in% l_locus[j]) %>% select(min.EHv312.a1) %>% pull() %>% as.integer() 
    aux_eh_a2 = val_data_superpopu %>% filter(locus %in% l_locus[j]) %>% select(max.EHv312.a2) %>% pull() %>% as.integer() 
    aux_eh_alleles_v2 = c(aux_eh_a1, aux_eh_a2)
    
    # Since we are keeping alleles for which sometimes an allele has EXP/NORMAL/PREMUT and not the other, we need to remove the same index
    index_na = which(is.na(aux_exp_alleles_v2))
    
    if (length(index_na) > 0){
      data_aux = xyTable(aux_exp_alleles_v2[-index_na], 
                         aux_eh_alleles_v2[-index_na])
    }else{
      data_aux = xyTable(aux_exp_alleles_v2[!is.na(aux_exp_alleles_v2)], 
                         aux_eh_alleles_v2[!is.na(aux_eh_alleles_v2)])
    }
    
    df_data_aux = data.frame(eh_alleles = data_aux$y,
                             exp_alleles = data_aux$x,
                             number_of_alleles = data_aux$number,
                             locus = rep(l_locus[j], length(data_aux$x)))
    
    
    # Concat info per locus
    df_data_with_freq_v2 = rbind(df_data_with_freq_v2,
                                 df_data_aux)
    
  }
  
  max_value = max(df_data_with_freq_v2$eh_alleles, 
                  df_data_with_freq_v2$exp_alleles,
                  na.rm = TRUE) + 5
  
  group.colors = rainbow(13)
  
  colour_locus = group.colors[i]
  
  superpopu = ggplot() +
    geom_point(data = df_data_with_freq_v2, 
               aes(x = exp_alleles, y = eh_alleles, size = number_of_alleles, color = factor(locus))) +
    xlim(5,max_value) +
    ylim(5,max_value) +
    geom_abline(method = "lm", formula = x ~ y, linetype = 2, colour = "gray") +  
    coord_equal() +
    scale_color_manual(values=group.colors) +
    labs(title = "", 
         y = "EH repeat sizes", 
         x = "PCR repeat sizes") + 
    #theme_light() +
    theme(legend.title = element_blank(),
          axis.title.x = element_blank()) +
    guides(size = FALSE) 
  
  superpopu_breakdown = ggplot() +
    geom_point(data = df_data_with_freq_v2, 
               aes(x = exp_alleles, y = eh_alleles, size = number_of_alleles, color = factor(locus))) +
    xlim(5,max_value) +
    ylim(5,max_value) +
    geom_abline(method = "lm", formula = x ~ y, linetype = 2, colour = "gray") +  
    coord_equal() +
    scale_color_manual(values=group.colors) +
    labs(title = "", 
         y = "EH repeat sizes", 
         x = "PCR repeat sizes") + 
    #theme_light() +
    theme(legend.title = element_blank(),
          axis.title.x = element_blank()) +
    guides(size = FALSE) +
    facet_wrap(locus~ .)
  
  output_superpopu = paste("Figure2/PCR_vs_EH_in", l_superpopus[i], sep = "_")
  output_superpopu_breakdown = paste(output_superpopu, "brokendown_by_locus", sep = "_")
  output_superpopu = paste(output_superpopu, ".png", sep = "")
  output_superpopu_breakdown = paste(output_superpopu_breakdown, ".png", sep = "")

  png(output_superpopu)
  print(superpopu)
  dev.off()
  
  png(output_superpopu_breakdown)
  print(superpopu_breakdown)
  dev.off()
}

