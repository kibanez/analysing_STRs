# Objective: create boxplots along with pvalue stats for case-control cases Arianna has date ()
Sys.info ()[c("nodename", "user")]
commandArgs ()
rm (list = ls ())
R.version.string ## "R version 3.6.1 (2019-07-05)"

# libraries
library(dplyr)
library(tidyverse)
library(ggplot2)
library(ggpubr)

# Set working directory
working_dir="~/Documents/STRs/ANALYSIS/cases_controls/Ari/"
setwd(working_dir)

# List all folders created by Ari so far
l_dir = list.files()

# we will go across all loci/folder
for (i in 1:length(l_dir)){
  # Let's read all case/control files
  l_files = list.files(l_dir[i])
  
  # create a new folder called `analysis`
  analysis_folder = paste(l_dir[i], "analysis", sep = "/")
  dir.create(analysis_folder)
  
  # read individual tables
  cases_motor_disorders = read.csv(paste(l_dir[i], "cases_motor_disorders_repeat_size.csv",sep = "/"), 
                                   sep = ",",
                                   stringsAsFactors = F,
                                   header = T)
  
  cases_neurodegen_motor_mitoch = read.csv(paste(l_dir[i], "cases_neurodegen_motor_mitoch_repeat_size.csv",sep = "/"), 
                                           sep = ",",
                                           stringsAsFactors = F,
                                           header = T)
  
  cases_neurodegenerative = read.csv(paste(l_dir[i], "cases_neurodegenerative_repeat_size.csv",sep = "/"), 
                                     sep = ",",
                                     stringsAsFactors = F,
                                     header = T)
  
  controls_adults_rd_cancer = read.csv(paste(l_dir[i], "controls_adults_rd_cancer_repeat_size.csv",sep = "/"), 
                                       sep = ",",
                                       stringsAsFactors = F,
                                       header = T)
  controls_adults_rd = read.csv(paste(l_dir[i], "controls_adults_rd_repeat_size.csv",sep = "/"), 
                                sep = ",",
                                stringsAsFactors = F,
                                header = T)
  
  
  # Merge all tables
  merged_table = rbind(cases_motor_disorders,
                       cases_neurodegen_motor_mitoch,
                       cases_neurodegenerative,
                       controls_adults_rd,
                       controls_adults_rd_cancer)
  
  
  # Generate boxplots with stats
  violin_plot = ggplot(merged_table, aes(x = group, y=repeat_size, fill = group)) +
    geom_violin() +
    xlab("Cases and controls datasets") + 
    ylab("Repeat sizes (repeat units)") + 
    theme(axis.text.x = element_text(angle = 90)) +
    ggtitle(l_dir[i]) 
  
  
  my_comparisons <- list( c("cases_motor_disorders_CNS", "controls_adult_rd"), c("cases_motor_disorders_CNS", "controls_adult_rd_cancer"),
                          c("cases_neurodeg", "controls_adult_rd"), c("cases_neurodeg", "controls_adult_rd_cancer"),
                          c("cases_neurodeg_motor_mito", "controls_adult_rd"), c("cases_neurodeg_motor_mito", "controls_adult_rd_cancer"))
  
  # Pairwise t-test between groups
  stat.test <- compare_means(repeat_size ~ group, 
                             data = merged_table,
                             paired = FALSE,
                             method = "wilcox.test",
                             p.adjust.method = "bonferroni",
                             alternative = "less") %>%
    mutate(y.position = c(51, 54, 57, 60, 63, 66, 69, 72, 75, 78))
  
  violin_plot_with_stat = ggviolin(merged_table, x = "group", y = "repeat_size", fill = "group") +
    stat_pvalue_manual(
      data = stat.test, label = "p.adj",
      #data = stat.test, label = "p.signif",
      xmin = "group1", xmax = "group2",
      y.position = "y.position") +
    theme(axis.text.x = element_text(angle = 90)) 
  
  violin_plot_with_stat_with_boxplot = ggviolin(merged_table, x = "group", y = "repeat_size", fill = "group",
                                                add = "boxplot", add.params = list(fill = "black")) +
    stat_pvalue_manual(
      data = stat.test, label = "p.adj",
      #data = stat.test, label = "p.signif",
      xmin = "group1", xmax = "group2",
      y.position = "y.position") +
    theme(axis.text.x = element_text(angle = 90)) 
  
  png(paste(l_dir[i], paste(paste("analysis/violin_plot_all_groups", l_dir[i], sep = "_"), ".png"),sep = "/"))
  print(violin_plot)
  dev.off()
  
  png(paste(l_dir[i], paste(paste("analysis/violin_plot_all_groups_with_stats", l_dir[i], sep = "_"), ".png"),sep = "/"))
  print(violin_plot_with_stat)
  dev.off()
  
  png(paste(l_dir[i], paste(paste("analysis/violin_plot_all_groups_with_stats_boxplot", l_dir[i], sep = "_"), ".png"),sep = "/"))
  print(violin_plot_with_stat_with_boxplot)
  dev.off()

  
  # Let's plot the same with ** or ns
  violin_plot_with_stat = ggviolin(merged_table, x = "group", y = "repeat_size", fill = "group") +
    stat_pvalue_manual(
      data = stat.test, label = "p.signif",
      xmin = "group1", xmax = "group2",
      y.position = "y.position") +
    theme(axis.text.x = element_text(angle = 90)) 
  
  violin_plot_with_stat_with_boxplot = ggviolin(merged_table, x = "group", y = "repeat_size", fill = "group",
                                                add = "boxplot", add.params = list(fill = "black")) +
    stat_pvalue_manual(
      data = stat.test, label = "p.signif",
      xmin = "group1", xmax = "group2",
      y.position = "y.position") +
    theme(axis.text.x = element_text(angle = 90)) 
  
  png(paste(l_dir[i], paste(paste("analysis/violin_plot_all_groups_with_ns", l_dir[i], sep = "_"), ".png"),sep = "/"))
  print(violin_plot_with_stat)
  dev.off()
  
  png(paste(l_dir[i], paste(paste("analysis/violin_plot_all_groups_with_ns_boxplot", l_dir[i], sep = "_"), ".png"),sep = "/"))
  print(violin_plot_with_stat_with_boxplot)
  dev.off()
  
}

