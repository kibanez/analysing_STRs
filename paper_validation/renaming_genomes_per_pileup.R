# Objective: rename the cases/genomes that present an expansion on GEL data with EHv2
# TP -> 60 alleles
# FP -> 10 alleles
date()
Sys.info()[c("nodename", "user")]
commandArgs()
rm(list = ls())
R.version.string # "R version 3.6.2 (2019-12-12)"

# Libraries
library(dplyr)

# Working directory
setwd("~/Documents/STRs/VALIDATION/QC_visual_inspection/pileups_june/")

# Load original or raw data
# Visual_QC_EHv2.5.5_ME.xlsx as tsv
l_pos_gel = read.csv("./list_GEL_EHv2_as_positive.txt", stringsAsFactors = F, header = F)
l_pos_gel = l_pos_gel$V1
length(l_pos_gel)
# 66

# Load correspondence between GEL Platekey and Paper ID
# From Table S7
df_ids = read.csv("./correspondence_gelID_paperID.tsv", stringsAsFactors = F, sep = "\t")
dim(df_ids)
# 635  3

