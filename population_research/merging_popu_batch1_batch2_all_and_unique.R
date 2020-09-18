# Objetive: merge population batch1 and batch2, in order to have a unique space
date ()
Sys.info ()[c("nodename", "user")]
commandArgs ()
rm (list = ls ())
R.version.string ## "R version 3.6.3 (2020-02-29)"

# libraries
library(dplyr); packageDescription ("dplyr", fields = "Version") #"0.8.5"
library(ggplot2); packageDescription ("ggplot2", fields = "Version") #"3.3.0"
library(reshape); packageDescription ("reshape", fields = "Version") #"0.8.8"
library(scatterplot3d); packageDescription("scatterplot3d", fields = "Version") # 0.3-41
library(ggpubr); packageDescription("ggpubr", fields = "Version") # 0.3.0

# Set environment
setwd("/Users/kibanez/Documents/STRs/ANALYSIS/population_research/MAIN_ANCESTRY/")

# For batch1, we are using Matthias' last work, with fine grained info
popu_batch1 = read.csv("./GEL_60k_germline_dataset_fine_grained_population_assignment20200224.csv",
                       stringsAsFactors = F,
                       header = T)
dim(popu_batch1)
# 59464  36

l_unrelated_batch1 = read.table("./60k_HWE_30k_random_unrelated_participants.txt", 
                                stringsAsFactors = F)
l_unrelated_batch1 = l_unrelated_batch1$V1
length(l_unrelated_batch1)
# 38344

# For batch2, we are using Thanos' work, and aggregated data in LK (latest release on RE, 3rd Sept)
popu_batch2 = read.csv("./batch2/aggV2_M30K_60K_1KGP3_ancestry_assignment_probs_R9_08062020.tsv",
                       stringsAsFactors = F,
                       header = T,
                       sep = " ")
dim(popu_batch2)
# 78388  33

l_unrelated_batch2 = read.table("./batch2/l_unrelated_55847_genomes_batch2.txt",
                                stringsAsFactors = F)
l_unrelated_batch2 = l_unrelated_batch2$V1
length(l_unrelated_batch2)
# 55847

# Check quality of batch2 comparing to batch1
l_intersected_genomes_b1_b2 = intersect(popu_batch1$ID, popu_batch2$plate_key)
length(l_intersected_genomes_b1_b2)
# 58,003

# are the predicted ancestry same?
# let's recode batch1
popu_batch1 = popu_batch1 %>%
  mutate(superpopu = case_when(best_guess_predicted_ancstry == "PJL" ~ "SAS",
                               best_guess_predicted_ancstry == "GBR" ~ "EUR",
                               best_guess_predicted_ancstry == "CEU" ~ "EUR",
                               best_guess_predicted_ancstry == "TSI" ~ "EUR",
                               best_guess_predicted_ancstry == "PUR" ~ "AMR",
                               best_guess_predicted_ancstry == "ACB" ~ "AFR",
                               best_guess_predicted_ancstry == "GIH" ~ "SAS",
                               best_guess_predicted_ancstry == "ASW" ~ "AFR",
                               best_guess_predicted_ancstry == "MXL" ~ "AMR",
                               best_guess_predicted_ancstry == "ESN" ~ "AFR",
                               best_guess_predicted_ancstry == "LWK" ~ "AFR",
                               best_guess_predicted_ancstry == "CHS" ~ "EAS",
                               best_guess_predicted_ancstry == "BEB" ~ "SAS",
                               best_guess_predicted_ancstry == "KHV" ~ "EAS",
                               best_guess_predicted_ancstry == "CLM" ~ "AMR",
                               best_guess_predicted_ancstry == "MSL" ~ "AFR",
                               best_guess_predicted_ancstry == "YRI" ~ "AFR",
                               best_guess_predicted_ancstry == "GWD" ~ "AFR",
                               best_guess_predicted_ancstry == "FIN" ~ "EUR",
                               best_guess_predicted_ancstry == "ITU" ~ "SAS",
                               best_guess_predicted_ancstry == "JPT" ~ "EAS",
                               best_guess_predicted_ancstry == "STU" ~ "SAS",
                               best_guess_predicted_ancstry == "CHB" ~ "EAS",
                               best_guess_predicted_ancstry == "CDX" ~ "EAS",
                               best_guess_predicted_ancstry == "PEL" ~ "AMR",
                               best_guess_predicted_ancstry == "IBS" ~ "EUR"))

# Do intersected genomes the same superpopu label in b1 and b2?
checking_intersected = left_join(popu_batch1 %>% filter(ID %in% l_intersected_genomes_b1_b2) %>% select(ID, best_guess_predicted_ancstry, self_reported, superpopu),
                                 popu_batch2 %>% filter(plate_key %in% l_intersected_genomes_b1_b2) %>% select(plate_key, ancestry0_8),
                                 by = c("ID" = "plate_key"))
dim(checking_intersected)
# 58003  5

checking_intersected = checking_intersected %>%
  mutate(are_equals = superpopu == ancestry0_8)
table(checking_intersected$are_equals)
#FALSE  TRUE 
#3979 54024

View(checking_intersected %>% filter(!are_equals))
# They are all FALSE == `unassigned` in batch2

l_not_consider_b2 = checking_intersected %>% filter(!are_equals) %>% select(ID) %>% pull() %>% unique()
length(l_not_consider_b2)
# 3979

# So, we can take ALL b1 + setidff(b2, l_not_consider_b2)
popu_merged = popu_batch1 %>% select(ID, superpopu)
aux_popu_batch2 = popu_batch2 %>% filter(!plate_key %in% l_not_consider_b2) %>% select(plate_key, ancestry0_8)
colnames(aux_popu_batch2) = c("ID", "superpopu")

popu_merged = rbind(popu_merged, aux_popu_batch2)
popu_merged = unique(popu_merged)
dim(popu_merged)
# 79849 2

# Union batch1 and batch2
# Since batch1 has been created from a random forest and is better curated, 