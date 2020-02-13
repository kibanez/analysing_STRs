# Objective: plot heatmap representing the corrected pvalue between different populations per forensics-loci
date ()
Sys.info ()[c("nodename", "user")]
commandArgs ()
rm (list = ls ())
R.version.string ## "R version 3.6.1 (2019-07-05)"

# libraries
library(dplyr)
library(ggplot2)
library(ggpubr)
library(reshape2)
library("RColorBrewer")

# Set working directory
working_dir="~/Documents/STRs/ANALYSIS/population_research/EH_3.1.2_research_October2019_55419_genomes/unrelated_probands_and_cancer/"
setwd(working_dir)

# Load main validation data table
main_data = read.csv("~/Documents/STRs/VALIDATION/EHv255_EHv312_validation_cohort_GEL_and_ILMN.tsv",
                     sep = "\t",
                     header = T,
                     stringsAsFactors = F)
l_loci = sort(unique(main_data$locus))

# remove PPP2R2B from there
l_loci = l_loci[-13]

# Create output folder
dir.create("heatmap")

# We need to do this by locus
# We want to produce an upper-triangle heatmap where, for each locus, the colour is associated with the corrected pvalue on
# how significant the repeat-size distribution in the locus is across different super-populations
for (i in 1:length(l_loci)){
  locus_name = l_loci[i]
  
  # AFR table
  afr_file = paste("AFR/cancer_and_RD/table_STR_repeat_size_each_row_allele_", locus_name, sep = "")
  afr_file = paste(afr_file, "_simplified_cancer_and_RD.tsv", sep = "")
  
  afr_table = read.csv(afr_file,
                       sep = "\t",
                       stringsAsFactors = F,
                       header = T)
  
  afr_table = afr_table %>%
    select(population, repeat_size)
  
  # AMR table
  amr_file = paste("AMR/cancer_and_RD/table_STR_repeat_size_each_row_allele_", locus_name, sep = "")
  amr_file = paste(amr_file, "_simplified_cancer_and_RD.tsv", sep = "")
  
  amr_table = read.csv(amr_file,
                       sep = "\t",
                       stringsAsFactors = F,
                       header = T)
  
  amr_table = amr_table %>%
    select(population, repeat_size)
  
  # ASI table
  asi_file = paste("ASI/cancer_and_RD/table_STR_repeat_size_each_row_allele_", locus_name, sep = "")
  asi_file = paste(asi_file, "_simplified_cancer_and_RD.tsv", sep = "")
  
  asi_table = read.csv(asi_file,
                       sep = "\t",
                       stringsAsFactors = F,
                       header = T)
  
  asi_table = asi_table %>%
    select(population, repeat_size)
  
  # EAS table
  eas_file = paste("EAS/cancer_and_RD/table_STR_repeat_size_each_row_allele_", locus_name, sep = "")
  eas_file = paste(eas_file, "_simplified_cancer_and_RD.tsv", sep = "")
  
  eas_table = read.csv(eas_file,
                       sep = "\t",
                       stringsAsFactors = F,
                       header = T)
  
  eas_table = eas_table %>%
    select(population, repeat_size)
  
  # EUR table
  eur_file = paste("EUR/cancer_and_RD/table_STR_repeat_size_each_row_allele_", locus_name, sep = "")
  eur_file = paste(eur_file, "_simplified_cancer_and_RD.tsv", sep = "")
  
  eur_table = read.csv(eur_file,
                       sep = "\t",
                       stringsAsFactors = F,
                       header = T)
  
  eur_table = eur_table %>%
    select(population, repeat_size)
  
  # Merged table
  merged_table = rbind(afr_table,
                       amr_table,
                       asi_table,
                       eas_table,
                       eur_table)
  
  # Pairwise t-test between groups
  stat.test <- compare_means(repeat_size ~ population, 
                             data = merged_table,
                             paired = FALSE,
                             method = "wilcox.test",
                             p.adjust.method = "bonferroni",
                             alternative = "great") %>%
    mutate(y.position = c(51, 53, 55, 57, 59, 61, 63, 65, 67, 69))
  
  # dataframe with padjusted values
  df_padj_locus = stat.test %>% select(group1, group2, p.adj) %>% as.data.frame()
  df_padj_locus$group1 = as.factor(df_padj_locus$group1)
  df_padj_locus$group2 = as.factor(df_padj_locus$group2)
  
  matrix_padj_locus = acast(df_padj_locus, group1 ~ group2, value.var = "p.adj")

  png(paste("heatmap", paste(locus_name, ".png", sep = ""), sep = "/"))
  heatmap.2(t(matrix_padj_locus),
            key = F,
            Rowv = F,
            trace = "none",
            dendrogram='none',
            main = locus_name)
  dev.off()
  
  #ggplot(data = matrix_padj_locus, 
  #       aes(x=group1, y=group2, fill=p.adj)) + 
  #  geom_tile(color = "white") +
  #  scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
  #                       midpoint = 0.5, limit = c(-1,1),  
  #                       name="Corrected\nWilcoxon test") +
  #  xlab("") + 
  #  ylab("") 
    
  
}
