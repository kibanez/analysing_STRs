# Function or R script, that given as input files PED and phenotypic files, it completes the PED file
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

# Reading arguments
options <- commandArgs(trailingOnly = T)

input_ped <- as.character(options[1])
genotype_ped <- as.character(options[2])
phenotype_file <- as.character(options[3])
output_ped <- as.character(options[4])

#Check whether `input_ped` and `phenotype_file` do exist
if (!file.exists(input_ped)){
  write("convert_VCF_into_PED R function: Original PED file does not exist. The original PED file is required", stderr())
}

if (!file.exists(phenotype_file)){
  write("convert_VCF_into_PED R function: Phenotypes (gender and affection status) for all genomes are missing", stderr())
}

# Load data
ped_data = read.csv(input_ped,
                    stringsAsFactors = F,
                    header = F,
                    sep = "\t")

pheno_data = read.csv(phenotype_file,
                      stringsAsFactors = F,
                      header = T,
                      sep = "\t")

genotype_data = read.csv(genotype_ped,
                         stringsAsFactors = F,
                         header = F,
                         sep = "\t")

# Enrich with `gender` and `affection status` each genome
ped_data = left_join(ped_data[,c(1:4)],
                     pheno_data %>% select(IID, sex, CaseControl),
                     by = c("V2" = "IID"))

# Merge `ped_data` with `genotypes`
ped_complete = cbind(ped_data,
                     genotype_data)

# Write into file
write.table(ped_complete,
            output_ped,
            quote = F, 
            row.names = F,
            col.names = F,
            sep = "\t")

