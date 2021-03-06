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

# Write each unrelated list of genomes for each super population

list_unrelated_EUR = popu_table %>%
  #filter(ID %in% list_unrelated_platekeys, self_reported %in% "European") %>%
  filter(ID %in% list_unrelated_platekeys, best_guess_predicted_ancstry %in% c("CEU", "FIN", "GBR", "IBS", "TSI")) %>%
  select(ID) %>%
  unique() %>%
  pull()
length(list_unrelated_EUR)
# 28153

list_unrelated_EAS = popu_table %>%
  filter(ID %in% list_unrelated_platekeys, best_guess_predicted_ancstry %in% c("KHV", "JPT", "CHS", "CHB")) %>%
  select(ID) %>%
  unique() %>%
  pull()
length(list_unrelated_EAS)
# 241

list_unrelated_AFR = popu_table %>%
  filter(ID %in% list_unrelated_platekeys, best_guess_predicted_ancstry %in% c("ACB","ASW", "GWD", "LWK", "MSL", "YRI", "ESN")) %>%
  select(ID) %>%
  unique() %>%
  pull()
length(list_unrelated_AFR)
# 1458

list_unrelated_AMR = popu_table %>%
  filter(ID %in% list_unrelated_platekeys, best_guess_predicted_ancstry %in% c("MXL", "PEL", "PUR", "CLM")) %>%
  select(ID) %>%
  unique() %>%
  pull()
length(list_unrelated_AMR)
# 965


list_unrelated_SAS = popu_table %>%
  filter(ID %in% list_unrelated_platekeys, best_guess_predicted_ancstry %in% c("PJL", "ITU", "GIH", "BEB", "STU")) %>%
  select(ID) %>%
  unique() %>%
  pull()
length(list_unrelated_SAS)
# 2897

# Write into files
write.table(list_unrelated_EUR, "list_EUR_28153_unrelated_genomes.txt", quote = F, col.names = F, row.names = F)
write.table(list_unrelated_EAS, "list_EAS_241_unrelated_genomes.txt", quote = F, col.names = F, row.names = F)
write.table(list_unrelated_AFR, "list_AFR_1458_unrelated_genomes.txt", quote = F, col.names = F, row.names = F)
write.table(list_unrelated_AMR, "list_AMR_965_unrelated_genomes.txt", quote = F, col.names = F, row.names = F)
write.table(list_unrelated_SAS, "list_SAS_2897_unrelated_genomes.txt", quote = F, col.names = F, row.names = F)

