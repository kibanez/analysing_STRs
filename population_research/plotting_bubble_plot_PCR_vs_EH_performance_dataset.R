# Objective: represent with bubble plots the correspondance between PCR and EH (after visual inspection) sizes
# split by super-population
date ()
Sys.info ()[c("nodename", "user")]
commandArgs ()
rm (list = ls ())
R.version.string ## "R version 3.6.3(2020-02-29)"

library(ggplot2); packageDescription ("ggplot2", fields = "Version") #"2.2.1"
library(grid); packageDescription ("grid", fields = "Version") #"3.3.2"
library(gridExtra); packageDescription ("gridExtra", fields = "Version") #"2.3"
library(reshape2); packageDescription ("reshape2", fields = "Version") #"1.4.3"
require(dplyr); packageDescription ("dplyr", fields = "Version") #"0.7.4"
library(lattice); packageDescription ("lattice", fields = "Version") #"0.20-35"
library(tidyverse); packageDescription ("tidyverse", fields = "Version") #"1.2.1"
library(RColorBrewer); packageDescription ("RColorBrewer", fields = "Version") #"1.1-2"
library(cowplot); packageDescription ("cowplot", fields = "Version") #"1.0.0"

# Set working environment
setwd("/Users/kibanez/Documents/STRs/PAPERS/POPULATION/figures/")

# Loading the golden or truthset or performance dataset
# Removing from here large expansions
# Only those < read-length
val_data = read.csv("~/Documents/STRs/data/PCR_EH/PCR_vs_EH_all_together_POPU_PAPER.tsv", sep = "\t", stringsAsFactors = F, header = T)
dim(val_data)
# 732  15

