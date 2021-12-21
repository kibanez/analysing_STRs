# Objective: analyse the insert-size across genomes within 100kGP
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
setwd("~/Documents/STRs/ANALYSIS/population_research/insert-size/")

# load data
agg_data = read.csv("./aggregate_gvcf_sample_stats_2021-11-19_16-50-11.tsv",
                    stringsAsFactors = F,
                    header = T,
                    sep = "\t")
dim(agg_data)
# 77646  67

# load the total number of genomes we are including in the population paper
l_unrel = read.table("./l_unrelated_55603_genomes_batch2.txt", stringsAsFactors = F)
l_unrel = l_unrel$V1
length(l_unrel)
# 55603

# load genomes sequenced at 125bp
l_125 = read.table("./list_genomes_125bp_100kGP.tsv", stringsAsFactors = F)
l_125 = l_125$V1
length(l_125)
# 15830

l_unrel_150 = setdiff(l_unrel, l_125)
length(l_unrel_150)
# 54437

agg_data = agg_data %>% 
  filter(platekey %in% l_unrel_150) %>%
  select(platekey, participant_id, samtools_insert_size_standard_deviation, samtools_insert_size_average) %>%
  unique()
dim(agg_data)
# 54288 4



