# Objective: for each locus/gene, plot altogether PCs vs repeat-sizes between different ancestries
date ()
Sys.info ()[c("nodename", "user")]
commandArgs ()
rm (list = ls ())
R.version.string ## "R version 3.6.2 (2019-12-12)"

# libraries
library(dplyr); packageDescription ("dplyr", fields = "Version") #"0.8.3"
library(ggplot2); packageDescription ("ggplot2", fields = "Version") #"3.2.1"
library("ggbeeswarm"); packageDescription ("ggbeeswarm", fields = "Version") #"0.6.0"
library("plot3D"); packageDescription ("plot3D", fields = "Version") #"1.3"



# Set environment
setwd("/Users/kibanez/Documents/STRs/ANALYSIS/population_research/EH_3.1.2_research_October2019_55419_genomes/unrelated_probands_and_cancer/")

# Load population data
popu_table_enriched = read.csv("../../population_info_enriched_59356_by_031019.tsv",
                               header = T,
                               sep = "\t",
                               stringsAsFactors = F)
dim(popu_table_enriched)
# 59356  20

# Take the list of all loci
# Load main validation data table
main_data = read.csv("~/Documents/STRs/VALIDATION/EHv255_EHv312_validation_cohort_GEL_and_ILMN.tsv",
                     sep = "\t",
                     header = T,
                     stringsAsFactors = F)
l_loci = sort(unique(main_data$locus))

# remove PPP2R2B from there
l_loci = l_loci[-13]

# Load merged repeat-sizes for each population
for(i in 1:length(l_loci)){
  locus_name = l_loci[i]
  
  print(locus_name)
  
  #AFR table
  afr_file = paste("AFR/cancer_and_RD/table_STR_repeat_size_each_row_allele_", locus_name, sep = "")
  afr_file = paste(afr_file, "_simplified_cancer_and_RD.tsv", sep = "")
  
  afr_table = read.csv(afr_file,
                       sep = "\t",
                       stringsAsFactors = F,
                       header = T)
  
  afr_table = afr_table %>%
    select(platekey, population, repeat_size)
  
  # AMR table
  amr_file = paste("AMR/cancer_and_RD/table_STR_repeat_size_each_row_allele_", locus_name, sep = "")
  amr_file = paste(amr_file, "_simplified_cancer_and_RD.tsv", sep = "")
  
  amr_table = read.csv(amr_file,
                       sep = "\t",
                       stringsAsFactors = F,
                       header = T)
  
  amr_table = amr_table %>%
    select(platekey, population, repeat_size)
  
  # ASI table
  asi_file = paste("ASI/cancer_and_RD/table_STR_repeat_size_each_row_allele_", locus_name, sep = "")
  asi_file = paste(asi_file, "_simplified_cancer_and_RD.tsv", sep = "")
  
  asi_table = read.csv(asi_file,
                       sep = "\t",
                       stringsAsFactors = F,
                       header = T)
  
  asi_table = asi_table %>%
    select(platekey, population, repeat_size)
  
  # EAS table
  eas_file = paste("EAS/cancer_and_RD/table_STR_repeat_size_each_row_allele_", locus_name, sep = "")
  eas_file = paste(eas_file, "_simplified_cancer_and_RD.tsv", sep = "")
  
  eas_table = read.csv(eas_file,
                       sep = "\t",
                       stringsAsFactors = F,
                       header = T)
  
  eas_table = eas_table %>%
    select(platekey, population, repeat_size)
  
  # EUR table
  eur_file = paste("EUR/cancer_and_RD/table_STR_repeat_size_each_row_allele_", locus_name, sep = "")
  eur_file = paste(eur_file, "_simplified_cancer_and_RD.tsv", sep = "")
  
  eur_table = read.csv(eur_file,
                       sep = "\t",
                       stringsAsFactors = F,
                       header = T)
  
  eur_table = eur_table %>%
    select(platekey, population, repeat_size)

  # Merged table
  merged_table = rbind(afr_table,
                       amr_table,
                       asi_table,
                       eas_table,
                       eur_table)
  
  # Enrich merged_table with PCs information
  merged_table = left_join(merged_table,
                           popu_table_enriched %>% select(platekey, pc1, pc2, pc3, pc4, pc5, pc6, pc7, pc8, pc9, pc10),
                           by = "platekey")
  
  pcs_vs_size_eur_afr = ggplot(data=merged_table %>% filter(population %in% c("EUR", "AFR")), 
         aes(x=pc2, y=repeat_size, colour = population)) +
    #geom_boxplot() +
    geom_violin() +
    geom_boxplot(width = 0.1) +
    #geom_boxplot(bins=300) +
    #geom_beeswarm(aes(color = population), grouponX=FALSE) +
    xlab("PC2") +
    ylab("Repeat-size") +
    guides(fill = FALSE)
  
  
  png(paste(paste("pc_vs_size/PCs_vs_size_EUR_vs_AFR_", locus_name, sep = ""), ".png", sep = ""))
  print(pcs_vs_size_eur_afr)
  dev.off()

  pcs_vs_size_eur_amr = ggplot(data=merged_table %>% filter(population %in% c("EUR", "AMR")), 
         aes(x=pc1, y=repeat_size, colour = population)) +
    #geom_boxplot() +
    geom_violin() +
    geom_boxplot(width = 0.1) +
    #geom_boxplot(bins=300) +
    #geom_beeswarm(aes(color = population), grouponX=FALSE) +
    xlab("PC1") +
    ylab("Repeat-size") +
    guides(fill = FALSE)
  
  png(paste(paste("pc_vs_size/PCs_vs_size_EUR_vs_AMR_", locus_name, sep = ""), ".png", sep = ""))
  print(pcs_vs_size_eur_amr)
  dev.off()

  pcs_vs_size_eur_eas = ggplot(data=merged_table %>% filter(population %in% c("EUR", "EAS")), 
         aes(x=pc1, y=repeat_size, colour = population)) +
    #geom_boxplot() +
    geom_violin() +
    geom_boxplot(width = 0.1) +
    #geom_boxplot(bins=300) +
    #geom_beeswarm(aes(color = population), grouponX=FALSE) +
    xlab("PC1") +
    ylab("Repeat-size") +
    guides(fill = FALSE)
  
  png(paste(paste("pc_vs_size/PCs_vs_size_EUR_vs_EAS_", locus_name, sep = ""), ".png", sep = ""))
  print(pcs_vs_size_eur_eas)
  dev.off()

  pcs_vs_size_eur_asi = ggplot(data=merged_table %>% filter(population %in% c("EUR", "ASI")), 
         aes(x=pc1, y=repeat_size, colour = population)) +
    #geom_boxplot() +
    geom_violin() +
    geom_boxplot(width = 0.1) +
    #geom_boxplot(bins=300) +
    #geom_beeswarm(aes(color = population), grouponX=FALSE) +
    xlab("PC2") +
    ylab("Repeat-size") +
    guides(fill = FALSE)
  
  png(paste(paste("pc_vs_size/PCs_vs_size_EUR_vs_ASI_", locus_name, sep = ""), ".png", sep = ""))
  print(pcs_vs_size_eur_asi)
  dev.off()

  pcs_vs_size_all = ggplot(data=merged_table, 
         aes(x=pc1, y=repeat_size, colour = population)) +
    #geom_boxplot() +
    geom_violin() +
    geom_boxplot(width = 0.2) +
    #geom_boxplot(bins=300) +
    #geom_beeswarm(aes(color = population), grouponX=FALSE) +
    xlab("PC1") +
    ylab("Repeat-size") +
    guides(fill = FALSE)
  
  png(paste(paste("pc_vs_size/PCs_vs_size_all_popus_", locus_name, sep = ""), ".png", sep = ""))
  print(pcs_vs_size_all)
  dev.off()
  
 # 3D
  #scatter3D(merged_table$pc1, 
  #          merged_table$pc2, 
  #          merged_table$repeat_size, 
  #          colvar = merged_table$population,
  #          phi = 0, bty ="g")
  
  
}
