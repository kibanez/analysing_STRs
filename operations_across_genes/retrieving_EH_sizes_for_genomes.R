# Objective: given a list of platekeys, we want to take the repeat-sizes for genes
date ()
Sys.info ()[c("nodename", "user")]
commandArgs ()
rm (list = ls ())
R.version.string ## "R version 3.3.2 (2016-10-31)"

# libraries
library(dplyr); packageDescription ("dplyr", fields = "Version") #"0.8.5"
library(tidyverse); packageDescription ("tidyverse", fields = "Version") # "1.2.1
library(ggplot2); packageDescription ("ggplot2", fields = "Version") #"3.3.0"

# set the working directory
setwd("~/Downloads/")

# Load the data
df_mito = read.csv("~/Downloads/table_mito_genomes.tsv",
                   stringsAsFactors = F,
                   header = T,
                   sep = "\t")
dim(df_mito)
# 345 2

# load merged august batch 2020
merged_table = read.csv("~/Documents/STRs/data/research/batch_august2020/output_EHv3.2.2_vcfs/merged/merged_93446_genomes_EHv322_batch_august2020.tsv",
                        stringsAsFactors = F,
                        header = T,
                        sep = "\t")
dim(merged_table)
# 27238  12

# Can you see if they are in your data for the population paper, and pull out a the repeat-sizes for all the 13 RE loci for all genomes? 
# A table with rows= genomes and columns =  26 (13 loci x 2 alleles)
list_genes = c("AR", "ATN1", "ATXN1", "ATXN2", "ATXN3", "ATXN7", "C9ORF72", "CACNA1A", "DMPK", "FMR1", "FXN","HTT", "TBP")
merged_table = merged_table %>%
  filter(gene %in% list_genes)



