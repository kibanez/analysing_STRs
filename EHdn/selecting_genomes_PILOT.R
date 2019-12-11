# Objective: select all PILOT genomes, but take their BAM file from frozen pilot dataset
date ()
Sys.info ()[c("nodename", "user")]
commandArgs ()
rm (list = ls ())
R.version.string ## "R version 3.6.1 (2019-07-05)"

# libraries
library(dplyr)
library(ggplot2)
library(reshape)

setwd("/Users/kibanez/Documents/STRs/ANALYSIS/EHdn/Pilot/")

pilot_data = read.csv("~/Documents/STRs/clinical_data/pilot_clinical_data/pilot_cohort_clinical_data_4833_genomes_removingPanels_280919.tsv",
                      header = T,
                      stringsAsFactors = F,
                      sep = "\t")
dim(pilot_data)
# 4974  10

paths_data= read.csv("~/Documents/STRs/clinical_data/pilot_clinical_data/b37_genomes_2019-12-11_11-07-20.tsv",
                     header = T,
                     stringsAsFactors = F, 
                     sep = "\t")
dim(paths_data)
# 4828  5

# How many genomes we do have in b37 absolute paths and genomes?
length(unique(paths_data$plate_key))
# 4828

# How many genomes we do have in frozen pilto data?
length(unique(pilot_data$plateKey))
# 4833

# Let's see at least that all 4824 are within 4833
length(intersect(unique(pilot_data$plateKey), unique(paths_data$plate_key)))
# 4828

df_pilot = paths_data %>% select(plate_key, path)
dim(df_pilot)
# 4828  2

for (i in 1:length(df_pilot$plate_key)){
  abs_path = paste(df_pilot$path[i], df_pilot$plate_key[i], sep = "/Assembly/")
  abs_path= paste(abs_path, "bam", sep = ".")
  df_pilot$path[i] = abs_path
}

# Write into final table separated by comma
write.table(df_pilot,
            "./EHdn-v0.8.6/input_EHdn-v0.8.6/list_4828_pilot_genomes_iSAAC.csv",
            quote = F,
            row.names = F,
            col.names = F,
            sep = ",")
