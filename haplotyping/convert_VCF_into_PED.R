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
phenotype_file <- as.character(options[2])
output_ped <- as.character(options[3])

#Check whether `input_ped` and `phenotype_file` do exist
if (!file.exists(input_ped)){
  write("convert_VCF_into_PED R function: Original PED file does not exist. The original PED file is required", stderr())
}

if (!file.exists(phenotype_file)){
  write("convert_VCF_into_PED R function: Phenotypes (gender and affection status) for all genomes are missing", stderr())
}


