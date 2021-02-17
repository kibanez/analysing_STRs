# Objective: create a barplot with the total number of unrel genomes in all datasets
date ()
Sys.info ()[c("nodename", "user")]
commandArgs ()
rm (list = ls ())
R.version.string ## "R version 3.6.1 (2019-07-05)"

# libraries
library(dplyr); packageDescription ("dplyr", fields = "Version") #"0.8.3"
library(ggplot2); packageDescription ("ggplot2", fields = "Version") #"3.2.1"
library(reshape); packageDescription ("reshape", fields = "Version") #"0.8.8"
library(scatterplot3d); packageDescription("scatterplot3d", fields = "Version") # 0.3-41
library(ggpubr); packageDescription("ggpubr", fields = "Version") # 0.2.3
library(tidyverse)

# Set working dir
setwd("~/Documents/STRs/ANALYSIS/population_research/")

# load raw data re total number of unrel genomes
total_unrel_genomes = read.csv("./number_genomes_per_superpopu_100k_1k_TOPMed.tsv",
                               stringsAsFactors = F, 
                               header = T,
                               sep = "\t")
dim(total_unrel_genomes)
# 15  3

barplot_all = ggplot(total_unrel_genomes, aes(x = population, y = number_genomes, fill = cohort)) +
  geom_bar(stat = "identity", position=position_dodge()) 

png("barplot_number_genomes_per_cohort_1kGP3_100kGP_TOPMed.png")
print(barplot_all)
dev.off()
