# Objective: take advantage from the work GEL bioinfo research group has done estimating ancestry for ~59,356 genomes
# and analyse and study the STR distribution across them
date ()
Sys.info ()[c("nodename", "user")]
commandArgs ()
rm (list = ls ())
R.version.string ## "R version 3.6.1 (2019-07-05)"

# libraries
library(dplyr)
library(ggplot2)
library(reshape)

setwd("/Users/kibanez/Documents/STRs/ANALYSIS/population_research/")


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
