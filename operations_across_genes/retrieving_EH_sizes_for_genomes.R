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


