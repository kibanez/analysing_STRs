# Objective: plot heatmap representing the corrected pvalue between different populations per disease-locus
# we will focus our disease-locus to those we have validated (13 as Feb'2020)
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
for (i in 1:length(l_loci)){
  
  
}