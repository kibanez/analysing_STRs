# Objective: analyse repeat-size distributions across different ancestries/super-populations and compare between all of them
date ()
Sys.info ()[c("nodename", "user")]
commandArgs ()
rm (list = ls ())
R.version.string ## "R version 3.6.1 (2019-07-05)"

# libraries
library(dplyr)
library(ggplot2)

# Set working directory
working_dir="~/Documents/STRs/ANALYSIS/population_research/EH_3.1.2_research_October2019_55419_genomes/unrelated_probands_and_cancer/"
setwd(working_dir)

# Load main data table
main_data = read.csv("./AFR/merged/merged_population_genomes_unrelated_probands_and_cancer_1136_avg_EHv3.1.2_AFR.tsv",
                     sep = "\t",
                     header = T,
                     stringsAsFactors = F)
l_loci = sort(unique(main_data$gene))

# We need to do this by locus
# Plot violing-plots as well

for (i in 1:length(l_loci)){

  # AFR table
  afr_file = paste("AFR/cancer_and_RD/table_STR_repeat_size_each_row_allele_", l_loci[i], sep = "")
  afr_file = paste(afr_file, "_simplified_cancer_and_RD.tsv", sep = "")
  
  afr_table = read.csv(afr_file,
                       sep = "\t",
                       stringsAsFactors = F,
                       header = T)
  
  afr_table = afr_table %>%
    select(population, repeat_size)
  
  # AMR table
  amr_file = paste("AMR/cancer_and_RD/table_STR_repeat_size_each_row_allele_", l_loci[i], sep = "")
  amr_file = paste(amr_file, "_simplified_cancer_and_RD.tsv", sep = "")
  
  amr_table = read.csv(amr_file,
                       sep = "\t",
                       stringsAsFactors = F,
                       header = T)
  
  amr_table = amr_table %>%
    select(population, repeat_size)
  
  # ASI table
  asi_file = paste("ASI/cancer_and_RD/table_STR_repeat_size_each_row_allele_", l_loci[i], sep = "")
  asi_file = paste(asi_file, "_simplified_cancer_and_RD.tsv", sep = "")
  
  asi_table = read.csv(asi_file,
                       sep = "\t",
                       stringsAsFactors = F,
                       header = T)
  
  asi_table = asi_table %>%
    select(population, repeat_size)
  
  # EAS table
  eas_file = paste("EAS/cancer_and_RD/table_STR_repeat_size_each_row_allele_", l_loci[i], sep = "")
  eas_file = paste(eas_file, "_simplified_cancer_and_RD.tsv", sep = "")
  
  eas_table = read.csv(eas_file,
                       sep = "\t",
                       stringsAsFactors = F,
                       header = T)
  
  eas_table = eas_table %>%
    select(population, repeat_size)
  
  # EUR table
  eur_file = paste("EUR/cancer_and_RD/table_STR_repeat_size_each_row_allele_", l_loci[i], sep = "")
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


  violin_plot = ggplot(merged_table, aes(x = population, y=repeat_size, fill = population)) +
    geom_violin() +
    xlab("Population") + 
    ylab("Repeat sizes (repeat units)") + 
    ggtitle(l_loci[i]) 
  
  
  my_comparisons <- list( c("AFR", "AMR"), c("AFR", "ASI"), c("AFR", "EAS"), c("AFR", "EUR"),
                          c("AMR", "ASI"), c("AMR", "EAS"), c("AMR", "EUR"))
  
  # Pairwise t-test between groups
  stat.test <- compare_means(repeat_size ~ population, 
                             data = merged_table,
                             paired = FALSE,
                             method = "wilcox.test",
                             p.adjust.method = "bonferroni",
                             alternative = "less") %>%
    mutate(y.position = c(51, 53, 55, 57, 59, 61, 63, 65, 67, 69))
  
  
  ggviolin(merged_table, x = "population", y = "repeat_size", fill = "population",
           add = "boxplot", add.params = list(fill = "black"))+
    stat_pvalue_manual(
      data = stat.test, label = "p.adj",
      xmin = "group1", xmax = "group2",
      y.position = "y.position") 
    #stat_compare_means(comparisons = my_comparisons, method = "wilcox.test", alternative = "less") + # Add significance levels
    #stat_compare_means(label.y = 50)                                       # Add global the p-value 
  
  
   png(paste(paste("analysis/violin_plot_all_ancestries", l_loci[i], sep = "_"), ".png"))
   print(violin_plot)
   dev.off()

}