# Objective: take advantage from the work GEL bioinfo research group has done estimating ancestry for ~59,356 genomes
# and analyse and study the STR distribution across them
date ()
Sys.info ()[c("nodename", "user")]
commandArgs ()
rm (list = ls ())
R.version.string ## "R version 3.6.1 (2019-07-05)"

# libraries
library(dplyr); packageDescription ("dplyr", fields = "Version") #"0.8.3"
library(ggplot2); packageDescription ("ggplot2", fields = "Version") #"3.2.1"
library(reshape); packageDescription ("reshape", fields = "Version") #"0.8.8"

setwd("/Users/kibanez/Documents/STRs/ANALYSIS/population_research/")

popu_table_enriched = read.csv("./population_info_enriched_59356_by_031019.tsv",
                               header = T,
                               sep = "\t",
                               stringsAsFactors = F)
dim(popu_table_enriched)
# 59356  20


ggplot(data=dat, aes(x=peddy_pc2, y=peddy_pc1)) +
  geom_hex(bins=100)

ggplot(data=dat, aes(x=peddy_pc2, y=peddy_pc1, colour=peddy_ancestry_pred)) +
  geom_hex(bins=100)

ggplot(data=dat, aes(x=peddy_pc3, y=peddy_pc1, colour=peddy_ancestry_pred)) +
  geom_hex(bins=100)

ggplot(data=dat, aes(x=peddy_pc4, y=peddy_pc1, colour=peddy_ancestry_pred)) +
  geom_hex(bins=100)

ggplot(data=dat, aes(x=peddy_pc3, y=peddy_pc2, colour=peddy_ancestry_pred)) +
  geom_hex(bins=100)

ggplot(data=dat, aes(x=peddy_pc4, y=peddy_pc2, colour=peddy_ancestry_pred)) +
  geom_hex(bins=100)

ggplot(data=dat, aes(x=peddy_pc3, y=peddy_pc4, colour=peddy_ancestry_pred)) +
  geom_hex(bins=100)

ggplot(data=dat, aes(x=1, y=peddy_ancestry_prob)) +
  geom_boxplot()

ggplot(data=dat, aes(x=peddy_ancestry_pred, y=peddy_ancestry_prob)) +
  geom_boxplot()
