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
working_dir="~/Documents/STRs/ANALYSIS/population_research/HipSTR/output_HipSTR/HipSTR_output_58971_vcfs/"
setwd(working_dir)

# Load main validation data table
main_data = read.csv("./EUR/merged/merged_forensics_loci_46883_EUR_HipSTRv0.6.2.tsv",
                     sep = "\t",
                     header = T,
                     stringsAsFactors = F)
l_loci = sort(unique(main_data$gene))

# Create output folder
dir.create("heatmap_forensics")

# We need to do this by locus
# We want to produce an upper-triangle heatmap where, for each locus, the colour is associated with the corrected pvalue on
# how significant the repeat-size distribution in the locus is across different super-populations
for (i in 1:length(l_loci)){
  locus_name = l_loci[i]
  
  # AFR table
  afr_table = read.csv("./AFR/merged/merged_forensics_loci_1777_AFR_HipSTRv0.6.2.tsv",
                       sep = "\t",
                       stringsAsFactors = F,
                       header = T)
  
  afr_table = afr_table %>%
    filter(gene %in% l_loci[i]) %>%
    select(repeat.size) %>%
    pull() %>%
    as.data.frame()
  colnames(afr_table) = "repeat.size"
  
  afr_table$population = rep("AFR", length(afr_table$repeat.size))
  
  # AMR table
  amr_table = read.csv("./AMR/merged/merged_forensics_loci_797_AMR_HipSTRv0.6.2.tsv",
                       sep = "\t",
                       stringsAsFactors = F,
                       header = T)
  
  amr_table = amr_table %>%
    filter(gene %in% l_loci[i]) %>%
    select(repeat.size) %>%
    pull() %>%
    as.data.frame()
  colnames(amr_table) = "repeat.size"
  
  amr_table$population = rep("AMR", length(amr_table$repeat.size))
  
  # ASI table
  asi_table = read.csv("./ASI/merged/merged_forensics_loci_5947_ASI_HipSTRv0.6.2.tsv",
                       sep = "\t",
                       stringsAsFactors = F,
                       header = T)
  
  asi_table = asi_table %>%
    filter(gene %in% l_loci[i]) %>%
    select(repeat.size) %>%
    pull() %>%
    as.data.frame()
  colnames(asi_table) = "repeat.size"
  
  asi_table$population = rep("ASI", length(asi_table$repeat.size))
  
  # EAS table
  eas_table = read.csv("./EAS/merged/merged_forensics_loci_400_EAS_HipSTRv0.6.2.tsv",
                       sep = "\t",
                       stringsAsFactors = F,
                       header = T)
  
  eas_table = eas_table %>%
    filter(gene %in% l_loci[i]) %>%
    select(repeat.size) %>%
    pull() %>%
    as.data.frame()
  colnames(eas_table) = "repeat.size"
  
  eas_table$population = rep("EAS", length(eas_table$repeat.size))
  
  
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
