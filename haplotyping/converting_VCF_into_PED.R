# Objective: convert multi-sample VCF file into PED file for further analysis using HaploView
date ()
Sys.info ()[c("nodename", "user")]
commandArgs ()
rm (list = ls ())
R.version.string ## "R version 3.6.1 (2019-07-05)"

# libraries
library(dplyr); packageDescription ("dplyr", fields = "Version") #"0.8.3"
library(ggplot2); packageDescription ("ggplot2", fields = "Version") #"3.2.1"
library(reshape); packageDescription ("reshape", fields = "Version") #"0.8.8"
library(scatterplot3d); packageDescription("scatterplot3d", fields = "Version") # 0.3-41
library(ggpubr); packageDescription("ggpubr", fields = "Version") # 0.2.3
library(tidyverse)

# Set environment
setwd("~/Documents/STRs/ANALYSIS/haplotyping/AR/HaploView/females/")

# Load data
ped_data = read.csv("./chrX_67495316-67595385.ped",
                    stringsAsFactors = F,
                    header = F,
                    sep = "\t")
dim(ped_data)
# 36  135

pheno_data = read.csv("./femalePhenotypeFile.txt",
                      stringsAsFactors = F,
                      header = T,
                      sep = " ")
dim(pheno_data)
# 38  3