# Objective: from the recoded and formatted tables from google drive (https://docs.google.com/spreadsheets/d/1cuh2rsDkQP3YEHjX6ogLWlgxO3sd0Jfs4zCVdizBQEc/edit#gid=1383787997)
# The aim here is to dedup all info we do have from 3 independent tables and merge them all
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

# Load 3 tables
gel_table = read.csv("./googleDrive_GEL_validation.tsv",
                     sep = "\t",
                     stringsAsFactors = F, 
                     header = T)
dim(gel_table)
# 635 11

ari_table = read.csv("./googleDrive_Arianna_NHNN_validation.tsv",
                     sep = "\t",
                     stringsAsFactors = F, 
                     header = T)
dim(ari_table)
# 144 11 

james_table = read.csv("./googleDrive_James_NHNN_validation.tsv",
                      sep = "\t",
                      stringsAsFactors = F, 
                      header = T)
dim(james_table)
# 48  11

merge_all = rbind(gel_table,
                  ari_table,
                  james_table)
dim(merge_all)
# 827  11

merge_all = unique(merge_all)
dim(merge_all)
# 820  11

# which are duplicates?
intersect(gel_table$LP_number, ari_table$LP_number)
# "LP3000999-DNA_C06" "LP3001031-DNA_H09" "LP3000329-DNA_E12" "LP3001101-DNA_H06" "LP3000118-DNA_E03" "LP3000124-DNA_D08"
intersect(gel_table$LP_number, james_table$LP_number)
# "LP3000595-DNA_E05" "LP3000469-DNA_C05" "LP3000474-DNA_E01"

# QC check -- they all have the same PCR< EHv2 and EHv3, estimations
write.table(merge_all, "googleDrive_all_merged_dedup_table.tsv", sep = "\t", quote = F, col.names = T, row.names = F)




