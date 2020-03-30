# Objective: enrich PCR_vs_EH estimation tables with population info
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

# Load Pilot programme with population info
pilot_popu = read.csv("~/Documents/STRs/ANALYSIS/population_research/PILOT_ANCESTRY/FINE_GRAINED_RF_classifications_incl_superPOP_prediction_final20191216.csv",
                      sep = ",",
                      stringsAsFactors = F,
                      header = T)
dim(pilot_popu)
# 4821  44

# Load Main programme with population info
main_popu = read.csv("~/Documents/STRs/ANALYSIS/population_research/MAIN_ANCESTRY/GEL_60k_germline_dataset_fine_grained_population_assignment20200224.csv",
                     sep = ",",
                     stringsAsFactors = F,
                     header = T)
dim(main_popu)
# 59464  36


# 1. GEL validation golden table
# I've copied the table as it is now (without super popu and subpopu info) from drive
# https://docs.google.com/spreadsheets/d/1cuh2rsDkQP3YEHjX6ogLWlgxO3sd0Jfs4zCVdizBQEc/edit#gid=0
gel_table = read.csv("./GEL_data_from_drive.tsv",
                     sep = "\t",
                     stringsAsFactors = F,
                     header = T)
dim(gel_table)
# 635  11

# Objective: enrich `Super.population` and `Sub,population` fields/columns (mix of Pilot and Main genomes)
gel_table = left_join(gel_table,
                      pilot_popu %>% select(ID, bestGUESS_sub_pop, bestGUESS_super_pop),
                      by = c("LP_number" = "ID"))
dim(gel_table)
# 635  13

# now with the main popu table
gel_table = left_join(gel_table,
                      main_popu %>% select(ID, best_guess_predicted_ancstry, self_reported),
                      by = c("LP_number" = "ID"))
dim(gel_table)
# 635  15

# Merge PILOT and MAIN columns into one
gel_table = gel_table %>%
  mutate(Super.population = coalesce(bestGUESS_super_pop, self_reported),
         Sub.population = coalesce(bestGUESS_sub_pop, best_guess_predicted_ancstry)) 

# Remove the extra final 4 columns
gel_table = gel_table[,-c(12:15)]
dim(gel_table)
# 635  11


# write into a file the final table
write.table(gel_table, 
            "./GEL_data_enriched_with_popu_to_put_in_drive.tsv",
            quote = F, 
            sep = "\t",
            row.names = F,
            col.names = T)


