# Objective: Select those family members (all within the same family) that have been recruited under `neuro` 
# Pilot and Main tables
# `disease_group` = neuro + `disease_subgroup` = Mitochondrial
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

# Set environment
setwd("/Users/kibanez/Documents/STRs/ANALYSIS/cases_controls/")

# Load PILOT clinical data
pilot_data = read.csv("~/Documents/STRs/clinical_data/pilot_clinical_data/pilot_cohort_clinical_data_4833_genomes_removingPanels_280919.tsv",
                      sep = "\t",
                      header = T,
                      stringsAsFactors = F)
dim(pilot_data)
# 4974  10

# Load V7 RE MAIN clinical data
main_data = read.csv("~/Documents/STRs/clinical_data/clinical_data/rd_genomes_all_data_250919.tsv",
                     sep = "\t",
                     header =  T,
                     stringsAsFactors = F)
dim(main_data)
# 1056568  26

# Load the table which translates different disease subgroups
translator_table = read.csv("~/Documents/STRs/clinical_data/clinical_research_cohort/phenotyping_v140_2019-09-13_15-26-02.tsv",
                            header = T,
                            sep = "\t",
                            stringsAsFactors = F)
dim(translator_table)
# 2632  3

