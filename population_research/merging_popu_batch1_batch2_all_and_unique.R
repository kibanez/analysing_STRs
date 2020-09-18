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

