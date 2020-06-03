# Objective: fish those UNRELATED genomes, already being annotated with population, that present an expanded STR in any of the 13 locus of the STR clinical paper
date ()
Sys.info ()[c("nodename", "user")]
commandArgs ()
rm (list = ls ())
R.version.string ## "R version 3.6.3 (2020-02-29)"

# libraries
library(dplyr)

# defining working directory
setwd("~/Documents/STRs/ANALYSIS/population_research/PAPER/expanded_genomes_main_pilot/")
