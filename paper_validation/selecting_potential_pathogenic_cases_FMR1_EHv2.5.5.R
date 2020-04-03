# Objective: from EHv2.5.5 summer batch select genomes with the following selection criteria (and compare with Arianna's numbers):
# both from GRCH37 and GRCH38:
# proband only
# specific disease == intellectual disability
# Arianna's numbers are the following:
# 5783 probands recruited under ID
# repeat size > 55
# 131 probands
date ()
Sys.info ()[c("nodename", "user")]
commandArgs ()
rm (list = ls ())
R.version.string ## "R version 3.6.1 (2019-07-05)"

library(ggplot2); packageDescription ("ggplot2", fields = "Version") #"3.2.1"
library(grid); packageDescription ("grid", fields = "Version") #"3.6.2"
library(gridExtra); packageDescription ("gridExtra", fields = "Version") #"2.3"
library(reshape2); packageDescription ("reshape2", fields = "Version") #"1.4.3"
require(dplyr); packageDescription ("dplyr", fields = "Version") #"0.8.3"

# Set working dir
setwd("~/Documents/STRs/data/research/EH_2.5.5_research_August2019/EH_output_v2.5.5_research_August_2019/merged_genotypeUpdated/")

# Load data
merged_table = read.csv("merged_loci_86457_research_genomes_new_loci_EHv2.5.5_summer2019_removingListVcfFiles.tsv",
                        sep = "\t",
                        stringsAsFactors = F,
                        header = T)
dim(merged_table)
# 3983  11

# Recode b37 and b38 chr names

