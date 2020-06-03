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

# Load data
main_table = read.csv("../../MAIN_ANCESTRY/GEL_60k_germline_dataset_fine_grained_population_assignment20200224.csv",
                      stringsAsFactors = F,
                      header = T,
                      sep = ",")
dim(main_table)
# 59464  36

pilot_table = read.csv("../../PILOT_ANCESTRY/FINE_GRAINED_RF_classifications_incl_superPOP_prediction_final20191216.csv",
                       stringsAsFactors = F,
                       header = T,
                       sep = ",")
dim(pilot_table)
# 4821  44

