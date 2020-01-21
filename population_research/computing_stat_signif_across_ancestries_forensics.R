# Objective: analyse expansions across forensics loci 
date ()
Sys.info ()[c("nodename", "user")]
commandArgs ()
rm (list = ls ())
R.version.string ## "R version 3.6.1 (2019-07-05)"

# libraries
library(dplyr)
library(ggplot2)
library(ggpubr)

# Set working directory
working_dir="~/Documents/STRs/ANALYSIS/population_research/HipSTR/output_HipSTR/HipSTR_output_58971_vcfs/"
setwd(working_dir)

# Load main data table
main_data = read.csv("./AFR/merged/merged_forensics_loci_1777_AFR_HipSTRv0.6.2.tsv",
                     sep = "\t",
                     header = T,
                     stringsAsFactors = F)
l_loci = sort(unique(main_data$gene))

# Create output folder
dir.create("analysis/")

# We need to do this by locus
# Plot violing-plots as well
for (i in 1:length(l_loci)){
  
  locus_name = l_loci[i]

  # AFR table
  afr_table = read.csv("AFR/merged/merged_forensics_loci_1777_AFR_HipSTRv0.6.2.tsv",
                       sep = "\t",
                       stringsAsFactors = F,
                       header = T)
  
  afr_table = afr_table %>%
    filter(gene %in% locus_name) %>%
    select(repeat.size)
  
  afr_table$population = rep("AFR", length(afr_table$repeat.size))
  
  # AMR table
  amr_table = read.csv("AMR/merged/merged_forensics_loci_797_AMR_HipSTRv0.6.2.tsv",
                       sep = "\t",
                       stringsAsFactors = F,
                       header = T)
  
  amr_table = amr_table %>%
    filter(gene %in% locus_name) %>%
    select(repeat.size)
  
  amr_table$population = rep("AMR", length(amr_table$repeat.size))
  
  # ASI table
  asi_table = read.csv("ASI/merged/merged_forensics_loci_5947_ASI_HipSTRv0.6.2.tsv",
                       sep = "\t",
                       stringsAsFactors = F,
                       header = T)
  
  asi_table = asi_table %>%
    filter(gene %in% locus_name) %>%
    select(repeat.size)
  
  asi_table$population = rep("ASI", length(asi_table$repeat.size))
  
  # EAS table
  eas_table = read.csv("EAS/merged/merged_forensics_loci_400_EAS_HipSTRv0.6.2.tsv",
                       sep = "\t",
                       stringsAsFactors = F,
                       header = T)
  
  eas_table = eas_table %>%
    filter(gene %in% locus_name) %>%
    select(repeat.size)
  
  eas_table$population = rep("EAS", length(eas_table$repeat.size))
  
  # EUR table
  eur_table = read.csv("EUR/merged/merged_forensics_loci_46883_EUR_HipSTRv0.6.2.tsv",
                       sep = "\t",
                       stringsAsFactors = F,
                       header = T)
  
  eur_table = eur_table %>%
    filter(gene %in% locus_name) %>%
    select(repeat.size)
  
  eur_table$population = rep("EUR", length(eur_table$repeat.size))
  
  # Merged table
  merged_table = rbind(afr_table,
                       amr_table,
                       asi_table,
                       eas_table,
                       eur_table)
  
  # Both the Mann-Whitney and the Kolmogorov-Smirnov tests are nonparametric tests to compare two unpaired groups of data. 
  # Both compute P values that test the null hypothesis that the two groups have the same distribution. But they work very differently:
  
  # The MANN-WHITNEY test first ranks all the values from low to high, and then computes a P value that depends on the discrepancy
  # between the mean ranks of the two groups.
  
  # The KOLMOGOROV-SMIRNOV test compares the cumulative distribution of the two data sets, and computes a P value that depends on the 
  # largest discrepancy between distributions.
  
#  Here are some guidelines for choosing between the two tests:
#  - The KS test is sensitive to any differences in the two distributions. Substantial differences in shape, spread or median will result in a small P value. 
  # In contrast, the MW test is mostly sensitive to changes in the median.
  
#  - The MW test is used more often and is recognized by more people, so choose it if you have no idea which to choose.
  
#  - The MW test has been extended to handle tied values. The KS test does not handle ties so well. If your data are categorical, so has many ties, don't choose the KS test.

#  - Some fields of science tend to prefer the KS test over the MW test. It makes sense to follow the traditions of your field.


  violin_plot = ggplot(merged_table, aes(x = population, y=repeat.size, fill = population)) +
    geom_violin() +
    xlab("Population") + 
    ylab("Repeat sizes (repeat units)") + 
    ggtitle(locus_name) 
  
  
  my_comparisons <- list( c("AFR", "AMR"), c("AFR", "ASI"), c("AFR", "EAS"), c("AFR", "EUR"),
                          c("AMR", "ASI"), c("AMR", "EAS"), c("AMR", "EUR"))
  
  # Pairwise t-test between groups
  stat.test <- compare_means(repeat.size ~ population, 
                             data = merged_table,
                             paired = FALSE,
                             method = "wilcox.test",
                             p.adjust.method = "bonferroni") %>%
                             #alternative = "great") %>%
    mutate(y.position = c(51, 53, 55, 57, 59, 61, 63, 65, 67, 69))
  
  violin_plot_with_stat = ggviolin(merged_table, x = "population", y = "repeat.size", fill = "population") +
    stat_pvalue_manual(
      #data = stat.test, label = "p.adj",
      data = stat.test, label = "p.signif",
      xmin = "group1", xmax = "group2",
      y.position = "y.position") 
  
  violin_plot_with_stat_with_boxplot = ggviolin(merged_table, x = "population", y = "repeat.size", fill = "population",
                                                add = "boxplot", add.params = list(fill = "black")) +
    stat_pvalue_manual(
      #data = stat.test, label = "p.adj",
      data = stat.test, label = "p.signif",
      xmin = "group1", xmax = "group2",
      y.position = "y.position") 
  
   png(paste(paste("analysis/violin_plot_all_ancestries", l_loci[i], sep = "_"), ".png"))
   print(violin_plot)
   dev.off()

   png(paste(paste("analysis/violin_plot_all_ancestries_with_stats", l_loci[i], sep = "_"), ".png"))
   print(violin_plot_with_stat)
   dev.off()
   
   png(paste(paste("analysis/violin_plot_all_ancestries_with_stats_boxplot", l_loci[i], sep = "_"), ".png"))
   print(violin_plot_with_stat_with_boxplot)
   dev.off()
   
}

