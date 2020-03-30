# Objective: Arianna worked out all PCR info retrieved from NHNN database
# I need to format the table
date ()
Sys.info ()[c("nodename", "user")]
commandArgs ()
rm (list = ls ())
R.version.string ## "R version 3.6.1 (2019-07-05)"

library(ggplot2); packageDescription ("ggplot2", fields = "Version") #"2.2.1"
library(grid); packageDescription ("grid", fields = "Version") #"3.3.2"
library(gridExtra); packageDescription ("gridExtra", fields = "Version") #"2.3"
library(reshape2); packageDescription ("reshape2", fields = "Version") #"1.4.3"
require(dplyr); packageDescription ("dplyr", fields = "Version") #"0.7.4"

# Set working directory
setwd("~/Documents/STRs/VALIDATION/PCR_EH_estimations/")

# Load Ari's NHNN PCR table
ari_table = read.csv("~/Documents/STRs/VALIDATION/Arianna_Fishing/STR Repeats in WinPath 06.03.20_per_KRI.tsv",
                     sep = "\t",
                     stringsAsFactors = F,
                     header = T)
dim(ari_table)
# 298  10

