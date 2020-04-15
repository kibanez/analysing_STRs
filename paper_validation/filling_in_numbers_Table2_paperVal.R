# Objective: R script to fill out Table2 of the paper

# libraries
library(dplyr)
library(tidyverse); packageDescription ("tidyverse", fields = "Version") # "1.2.1

# Set directory
setwd("~/Documents/STRs/PAPERS/VALIDATION_PAPER/")

# Load main table
table_diseases = read.csv("table_diseases_enriched_popu_includingSkeletalMuscleChan.tsv",
                          stringsAsFactors = F, 
                          header = T,
                          sep = "\t")
dim(table_diseases)
# 11801  16

# Load pilot table
