# Objective: create bubble plots representing repeat sizes EH vs PCR, colouring by super-population
# PER LOCUS
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

# Load data
merged_table = read.csv("./googleDrive_all_merged_dedup_table.tsv",
                        sep = "\t",
                        stringsAsFactors = F, 
                        header = T)
dim(merged_table)
# 792  11

merged_table$min.PCR.a1 = as.numeric(merged_table$min.PCR.a1)
merged_table$maxPCR.a2 = as.numeric(merged_table$maxPCR.a2)
merged_table$min.EHv255.a1 = as.numeric(merged_table$min.EHv255.a1)
merged_table$max.EHv255.a2 = as.numeric(merged_table$max.EHv255.a2)
merged_table$min.EHv312.a1 = as.numeric(merged_table$min.EHv312.a1)
merged_table$max.EHv312.a2 = as.numeric(merged_table$max.EHv312.a2)

# Per locus, we want to plot a bubble plot represeting the repeat-sizes from EHv2/EHv3 vs PCR
# and colouring by diff colours superpopu 
l_loci = unique(merged_table$locus)

for (i in 1:length(l_loci)){
  aux_table = merged_table %>%
    filter(locus %in% l_loci[i])
  
  exp_alleles = c(aux_table$min.PCR.a1, aux_table$maxPCR.a2)
  ehv2_alleles = c(aux_table$min.EHv255.a1, aux_table$max.EHv255.a2)
  ehv3_alleles = c(aux_table$min.EHv312.a1, aux_table$max.EHv312.a2)
  superpopu = c(aux_table$Super.population, aux_table$Super.population)
  
  df_to_plot = data.frame(exp_alleles, ehv2_alleles, ehv3_alleles, superpopu)
  df_to_plot$superpopu = as.character(df_to_plot$superpopu)
  
  max_valuev2 = max(exp_alleles, ehv2_alleles, na.rm = T)
  max_valuev3 = max(exp_alleles, ehv3_alleles, na.rm = T)
  
  # Remove NA's, because we are not going to plot them
  df_to_plot = df_to_plot[complete.cases(df_to_plot),]
  
  # Jitter, to avoid overlapping superpopu points
  #set.seed(3)
  # Vary the marker size
  
  locus_bubble_v2 = ggplot(df_to_plot, 
                        aes(x = exp_alleles, y = ehv2_alleles, colour = superpopu, size=superpopu)) + 
      geom_point() + 
    scale_color_manual(values=c("red","green","blue","purple", "yellow")) +
    scale_size_manual(values=c(9,7,5,3,1)) +
      xlim(5,max_valuev2) + 
      ylim(5,max_valuev2) + 
      labs(title = "", 
           y = "Repeat sizes for each allele \n EHv2.5.5", 
           x = "Repeat sizes for each allele \n PCR") + 
      coord_equal() +
      theme(legend.title = element_blank()) + 
      guides(size = FALSE)
  
  locus_bubble_v3 = ggplot(df_to_plot, 
                           aes(x = exp_alleles, y = ehv3_alleles, colour = superpopu, size=superpopu)) + 
    geom_point() + 
    scale_color_manual(values=c("red","green","blue","purple", "yellow")) +
    scale_size_manual(values=c(9,7,5,3,1)) +
    xlim(5,max_valuev2) + 
    ylim(5,max_valuev2) + 
    labs(title = "", 
         y = "Repeat sizes for each allele \n EHv3.1.2", 
         x = "Repeat sizes for each allele \n PCR") + 
    coord_equal() +
    theme(legend.title = element_blank()) + 
    guides(size = FALSE)
  
  output_ehv2_locus = paste("./figures/", paste(l_loci[i], "PCR_EHv255_across_superpopup.png", sep = "_"), sep = "")
  output_ehv3_locus = paste("./figures/", paste(l_loci[i], "PCR_EHv312_across_superpopup.png", sep = "_"), sep = "")
  
  png(output_ehv2_locus)
  print(locus_bubble_v2)
  dev.off()
  
  png(output_ehv3_locus)
  print(locus_bubble_v3)
  dev.off()

}


locus_bubble_v3 = ggplot(df_to_plot, 
                         aes(x = exp_alleles, y = ehv3_alleles, colour = superpopu)) + 
  geom_point(size = 2, aes(fill = superpopu)) + 
  xlim(5,max_valuev3) + 
  ylim(5,max_valuev3) + 
  labs(title = "", 
       y = "Repeat sizes for each allele \n EHv3.1.2", 
       x = "Repeat sizes for each allele \n PCR") + 
  coord_equal() +
  theme(legend.title = element_blank()) + 
  guides(size = FALSE)





