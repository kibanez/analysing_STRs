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



