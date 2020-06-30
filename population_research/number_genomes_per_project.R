# Objective: combine 100K, 1K, and gnomAD projects and see how many genomes we are analysing
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
setwd("~/Documents/STRs/ANALYSIS/population_research/PAPER/")

# gnomAD
# Load data
gnomad_popu = read.csv("../gnomAD/number_genomes_per_superpopu.tsv",
                       stringsAsFactors = F,
                       header = F,
                       sep = "\t")
dim(gnomad_popu)
# 9  2

colnames(gnomad_popu) = c("superpopu", "number")

# Load data only CRAMs
gnomad_popu_crams = read.csv("../gnomAD/number_genomes_per_superpopu_with_CRAMs.tsv",
                       stringsAsFactors = F,
                       header = F,
                       sep = "\t")
dim(gnomad_popu_crams)
# 9  2
colnames(gnomad_popu_crams) = c("superpopu", "number")


png("../gnomAD/figures/barplot_ancestries_groups_raw_numbers.png")
ggplot(gnomad_popu, 
       aes(x = reorder(superpopu, - number), y = number)) + 
  geom_bar(stat = "identity", aes(fill = superpopu)) + 
  geom_text(aes(label=number), vjust=-0.5, size = 4, colour = "grey") +
  ylab("Number of genomes") + 
  xlab("Ancestry cohorts - gnomAD V3 -  29,074 total genomes") 
dev.off()

# with crams
png("../gnomAD/figures/barplot_ancestries_groups_raw_numbers_with_CRAMs.png")
ggplot(gnomad_popu_crams, 
       aes(x = reorder(superpopu, - number), y = number)) + 
  geom_bar(stat = "identity", aes(fill = superpopu)) + 
  geom_text(aes(label=number), vjust=-0.5, size = 4, colour = "grey") +
  ylab("Number of genomes") + 
  xlab("Ancestry cohorts - gnomAD V3 -  16,392 total genomes (whole CRAMs)") 
dev.off()

# 1K cohort
# load data
thousand_data = read.csv("../1kg/1000G_2504_high_coverage.sequence.index.tsv",
                             stringsAsFactors = F,
                             header = T,
                             sep = "\t")
dim(thousand_data)
# 2504 22

thousand_data = thousand_data %>%
  select(POPULATION, SAMPLE_NAME)
thousand_data = unique(thousand_data)
dim(thousand_data)
# 2504  2

# recode superpopu
thousand_data = thousand_data %>%
  mutate(superpopu = case_when(POPULATION == "PJL" ~ "SAS",
                               POPULATION == "GBR" ~ "EUR",
                               POPULATION == "CEU" ~ "EUR",
                               POPULATION == "TSI" ~ "EUR",
                               POPULATION == "PUR" ~ "AMR",
                               POPULATION == "ACB" ~ "AFR",
                               POPULATION == "GIH" ~ "SAS",
                               POPULATION == "ASW" ~ "AFR",
                               POPULATION == "MXL" ~ "AMR",
                               POPULATION == "ESN" ~ "AFR",
                               POPULATION == "LWK" ~ "AFR",
                               POPULATION == "CHS" ~ "EAS",
                               POPULATION == "BEB" ~ "SAS",
                               POPULATION == "KHV" ~ "EAS",
                               POPULATION == "CLM" ~ "AMR",
                               POPULATION == "MSL" ~ "AFR",
                               POPULATION == "YRI" ~ "AFR",
                               POPULATION == "GWD" ~ "AFR",
                               POPULATION == "FIN" ~ "EUR",
                               POPULATION == "ITU" ~ "SAS",
                               POPULATION == "JPT" ~ "EAS",
                               POPULATION == "STU" ~ "SAS",
                               POPULATION == "CHB" ~ "EAS",
                               POPULATION == "PEL" ~ "AMR",
                               POPULATION == "IBS" ~ "EUR"))

thousand_data = thousand_data %>%
  group_by(superpopu) %>%
  mutate(number = n()) %>%
  ungroup() %>%
  as.data.frame()

thousand_data_numbers = thousand_data %>%
  select(superpopu, number)
thousand_data_numbers = unique(thousand_data_numbers)

png("../1kg/figures/barplot_ancestries_groups_raw_numbers.png")
ggplot(thousand_data_numbers, 
       aes(x = reorder(superpopu, - number), y = number)) + 
  geom_bar(stat = "identity", aes(fill = superpopu)) + 
  geom_text(aes(label=number), vjust=-0.5, size = 4, colour = "grey") +
  ylab("Number of genomes") + 
  xlab("Ancestry cohorts - 1K Phase3 -  2,504 total genomes") 
dev.off()
