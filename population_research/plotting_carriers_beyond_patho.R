# Objective: plot population distribution for genomes presenting a carrier beyond full mutation or pathogenic cutoff
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

# Set working dir
setwd("~/Documents/STRs/ANALYSIS/population_research/100K/carrier_freq/")

# Load MAIN popu table we have so far
popu_table = read.csv("~/Documents/STRs/ANALYSIS/population_research/MAIN_ANCESTRY/GEL_60k_germline_dataset_fine_grained_population_assignment20200224.csv",
                      stringsAsFactors = F, 
                      sep = ",",
                      header = T)
dim(popu_table)
# 59464  36

# Load PILOT popu table 
pilot_popu_table = read.csv("~/Documents/STRs/ANALYSIS/population_research/PILOT_ANCESTRY/FINE_GRAINED_RF_classifications_incl_superPOP_prediction_final20191216.csv",
                            stringsAsFactors = F,
                            sep = ",",
                            header = T)
dim(pilot_popu_table)
# 4821  44 

# Load HTT EHv322 batch march 2020, after visual QC
htt_table = read.csv("./list_PIDs_for_HTT_pileup.tsv",
                    stringsAsFactors = F,
                    sep = "\t",
                    header = T)
dim(htt_table)
# 231  4

l_htt = htt_table %>%
  filter(larger.than.40.after.visual.QC. %in% "yes") %>%
  select(PLATEKEY) %>%
  unique() %>% 
  pull()
length(l_htt)
# 51

