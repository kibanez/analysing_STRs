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

popu_table = popu_table %>%
  mutate(merged_superpopu = case_when(best_guess_predicted_ancstry == "ACB" ~ "AFR",
                                      best_guess_predicted_ancstry == "ASW" ~ "AFR",
                                      best_guess_predicted_ancstry == "BEB" ~ "SAS",
                                      best_guess_predicted_ancstry == "CEU" ~ "EUR",
                                      best_guess_predicted_ancstry == "CHB" ~ "EAS",
                                      best_guess_predicted_ancstry == "CHS" ~ "EAS",
                                      best_guess_predicted_ancstry == "CLM" ~ "AMR",
                                      best_guess_predicted_ancstry == "ESN" ~ "AFR",
                                      best_guess_predicted_ancstry == "FIN" ~ "EUR",
                                      best_guess_predicted_ancstry == "GBR" ~ "EUR",
                                      best_guess_predicted_ancstry == "GIH" ~ "SAS",
                                      best_guess_predicted_ancstry == "GWD" ~ "AFR",
                                      best_guess_predicted_ancstry == "IBS" ~ "EUR",
                                      best_guess_predicted_ancstry == "ITU" ~ "SAS",
                                      best_guess_predicted_ancstry == "JPT" ~ "EAS",
                                      best_guess_predicted_ancstry == "KHV" ~ "AFR",
                                      best_guess_predicted_ancstry == "LWK" ~ "AFR",
                                      best_guess_predicted_ancstry == "MSL" ~ "AFR",
                                      best_guess_predicted_ancstry == "MXL" ~ "AMR",
                                      best_guess_predicted_ancstry == "PEL" ~ "AMR",
                                      best_guess_predicted_ancstry == "PJL" ~ "SAS",
                                      best_guess_predicted_ancstry == "PUR" ~ "AMR",
                                      best_guess_predicted_ancstry == "STU" ~ "SAS",                                      
                                      best_guess_predicted_ancstry == "TSI" ~ "EUR",
                                      best_guess_predicted_ancstry == "YRI" ~ "AFR"))

# Load PILOT popu table 
pilot_popu_table = read.csv("~/Documents/STRs/ANALYSIS/population_research/PILOT_ANCESTRY/FINE_GRAINED_RF_classifications_incl_superPOP_prediction_final20191216.csv",
                            stringsAsFactors = F,
                            sep = ",",
                            header = T)
dim(pilot_popu_table)
# 4821  44 

pilot_popu_table = pilot_popu_table %>%
  mutate(merged_superpopu_pilot = case_when(bestGUESS_sub_pop == "ACB" ~ "AFR",
                                            bestGUESS_sub_pop == "ASW" ~ "AFR",
                                            bestGUESS_sub_pop == "BEB" ~ "SAS",
                                            bestGUESS_sub_pop == "CEU" ~ "EUR",
                                            bestGUESS_sub_pop == "CHB" ~ "EAS",
                                            bestGUESS_sub_pop == "CHS" ~ "EAS",
                                            bestGUESS_sub_pop == "CLM" ~ "AMR",
                                            bestGUESS_sub_pop == "ESN" ~ "AFR",
                                            bestGUESS_sub_pop == "GBR" ~ "EUR",
                                            bestGUESS_sub_pop == "GIH" ~ "SAS",
                                            bestGUESS_sub_pop == "GWD" ~ "AFR",
                                            bestGUESS_sub_pop == "IBS" ~ "EUR",
                                            bestGUESS_sub_pop == "ITU" ~ "SAS",
                                            bestGUESS_sub_pop == "KHV" ~ "AFR",
                                            bestGUESS_sub_pop == "LWK" ~ "AFR",
                                            bestGUESS_sub_pop == "MSL" ~ "AFR",
                                            bestGUESS_sub_pop == "MXL" ~ "AMR",
                                            bestGUESS_sub_pop == "PJL" ~ "SAS",
                                            bestGUESS_sub_pop == "PUR" ~ "AMR",
                                            bestGUESS_sub_pop == "STU" ~ "SAS",
                                            bestGUESS_sub_pop == "TSI" ~ "EUR",
                                            bestGUESS_sub_pop == "YRI" ~ "AFR"))

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

popu_htt = popu_table %>%
  filter(ID %in% l_htt) %>%
  select(ID, merged_superpopu)
popu_htt = unique(popu_htt)
dim(popu_htt)
# 26  2


pilot_popu_htt = pilot_popu_table %>%
  filter(ID %in% l_htt) %>%
  select(ID, merged_superpopu_pilot)
pilot_popu_htt = unique(pilot_popu_htt)
dim(pilot_popu_htt)
# 1  2
colnames(pilot_popu_htt) = colnames(popu_htt)

merged_popu_htt = rbind(popu_htt,
                        pilot_popu_htt)
dim(merged_popu_htt)
# 27  2

