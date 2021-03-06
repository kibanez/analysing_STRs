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
setwd("~/Documents/STRs/ANALYSIS/haplotyping/AR/HaploView/2021/together_female_male/")

# Load data
ped_data = read.csv("./chrX_67495316-67595385.ped",
                    stringsAsFactors = F,
                    header = F,
                    sep = "\t")
dim(ped_data)
# 71 11

pheno_data = read.csv("./female_and_male_PhenotypeFile.txt",
                      stringsAsFactors = F,
                      header = T,
                      sep = " ")
dim(pheno_data)
# 75  4

# Enrich with `gender` and `affection status` each genome
ped_data = left_join(ped_data,
                     pheno_data %>% select(IID, sex, CaseControl),
                     by = c("V2" = "IID"))

# In this case specify 5th column as sex
ped_data$V5 = ped_data$sex
# Specify affection status
ped_data$V6 = ped_data$CaseControl


# Write final PED file into a file
# The first 6 columns are TAB separated
# The genotypes are space and TAB separated (<allele1 allele2>\t<allele1 allele2>)
write.table(ped_data %>% select(V1,V2,V3,V4,V5,V6),
            "./together_female_male_cc_header.ped",
            quote = F,
            row.names = F,
            col.names = F,
            sep = "\t")

# Merge `female_cc_AF3_header.ped` and `chrX_67495316-67595385_genotypes.ped`
only_genotypes = read.csv("./chrX_67495316-67595385_genotypes.ped",
                          sep = "\t",
                          stringsAsFactors = F,
                          header = F)
dim(only_genotypes)
# 71  5

merged_final_ped = cbind(ped_data %>% select(V1,V2,V3,V4,V5,V6),
                         only_genotypes)
write.table(merged_final_ped,
            "./chrX_67495316-67595385_complete.ped",
            quote = F,
            row.names = F,
            col.names = F,
            sep = "\t")
