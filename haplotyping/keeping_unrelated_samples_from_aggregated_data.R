# Objective: for haplotyping, any locus or region, we need to work with unrelated samples or genomes
# THe aim here is to take unrelated genomes from Loukas group's work on population (Main Programme)
date ()
Sys.info ()[c("nodename", "user")]
commandArgs ()
rm (list = ls ())
R.version.string ## "R version 3.6.1 (2019-07-05)"

library(ggplot2); packageDescription ("ggplot2", fields = "Version") #"3.2.1"
library(grid); packageDescription ("grid", fields = "Version") #"3.6.1"
library(gridExtra); packageDescription ("gridExtra", fields = "Version") #"2.3"
library(reshape2); packageDescription ("reshape2", fields = "Version") #"1.4.3"
require(dplyr); packageDescription ("dplyr", fields = "Version") #"0.8.3"

# Working directory
setwd("~/Documents/STRs/ANALYSIS/haplotyping/")

# load main data
popu_table = read.csv("~/Documents/STRs/ANALYSIS/population_research/MAIN_ANCESTRY/GEL_60k_germline_dataset_fine_grained_population_assignment20200224.csv",
                      stringsAsFactors = F, 
                      header= T)
dim(popu_table)
# 59464  36


# clinical data
clin_data = read.csv("~/Documents/STRs/clinical_data/clinical_data/rd_genomes_all_data_300320.tsv",
                     stringsAsFactors = F, 
                     header = T,
                     sep = "\t")
dim(clin_data)
# 1124633  31

popu_table = left_join(popu_table,
                       clin_data %>% select(platekey, rare_diseases_family_id, participant_type, programme, participant_stated_gender),
                       by = c("ID" = "platekey"))
popu_table = unique(popu_table)
dim(popu_table)
# 59464  40

# List of unrelated = probands + is.na(participant_type)
list_unrelated = unique(c(which(popu_table$participant_type %in% "Proband"),
                   which(is.na(popu_table$participant_type)),
                   which(popu_table$programme %in% "Cancer")))
list_unrelated_platekeys = popu_table$ID[list_unrelated]
length(list_unrelated_platekeys)
# 33714

# Write it into a file
write.table(list_unrelated_platekeys, "list_33714_unrelated_genomes.txt", quote = F, col.names = F, row.names = F)
