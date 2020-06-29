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

# Load data
gnomad_popu = read.csv("../gnomAD/number_genomes_per_superpopu.tsv",
                       stringsAsFactors = F,
                       header = F,
                       sep = "\t")
dim(gnomad_popu)
# 9  2

colnames(gnomad_popu) = c("superpopu", "number")

png("../gnomAD/figures/barplot_ancestries_groups_raw_numbers.png")
ggplot(gnomad_popu, 
       aes(x = reorder(superpopu, - number), y = number)) + 
  geom_bar(stat = "identity", aes(fill = superpopu)) + 
  geom_text(aes(label=number), vjust=-0.5, size = 4, colour = "grey") +
  ylab("Number of genomes") + 
  xlab("Ancestry cohorts - gnomAD V3 -  29,074 total genomes") 
dev.off()
