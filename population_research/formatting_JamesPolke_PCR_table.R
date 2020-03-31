# Objective: James' team worked through some extra genomes, estimating PCR repeat-size alleles
# I need to format the table and enrich with super and sub-population and gender information
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

# Load James Polke TSV file (worked from the original xls file)
james_table = read.csv("../JamesPolke/jamesPolke_non_NA_results.tsv",
                       sep = "\t",
                       stringsAsFactors = F,
                       header = T)
dim(james_table)
# 48  4

# First take the latest platekey corresponding to the PID
clin_data = read.csv("~/Documents/STRs/clinical_data/clinical_data/rd_genomes_all_data_300320.tsv",
                     sep = "\t",
                     stringsAsFactors = F,
                     header = T)
dim(clin_data)
# 1124633  31

# population table
popu_table = read.csv("~/Documents/STRs/ANALYSIS/population_research/MAIN_ANCESTRY/GEL_60k_germline_dataset_fine_grained_population_assignment20200224.csv",
                      sep = ",",
                      stringsAsFactors = F,
                      header = T)
dim(popu_table)
# 59464  36

# Corresponding platekey <-> PID
james_table = left_join(james_table,
                        clin_data %>% select(participant_id, plate_key),
                        by = c("PID" = "participant_id"))
james_table = unique(james_table)
dim(james_table)
# 48 5

# Write into a file platekey, locus, in order to retrieve EHv2/EHv3 estimations for them
to_retrieve = james_table %>% select(plate_key, locus)
write.table(to_retrieve, "../JamesPolke/platekey_locus_JamesPolke.csv", sep = ",", quote = F, row.names = F, col.names = F)

# Once retrieved EHv2/EHv3 info, let's take this
james_ehv2_ehv3 = read.csv("../JamesPolke/platekey_locus_JamesPolke_EHv2_EHv3.tsv",
                           sep = "\t",
                           stringsAsFactors = F, 
                           header = T)
dim(james_ehv2_ehv3)
# 48  6

james_all = left_join(james_table,
                      james_ehv2_ehv3,
                      by = c("plate_key" = "platekey", "locus" = "locus"))
dim(james_all)
# 48  9


# Let's enrich with population
james_all = left_join(james_all,
                      popu_table %>% select(ID, best_guess_predicted_ancstry, self_reported),
                      by = c("plate_key" = "ID"))
dim(james_all)
# 48  11


