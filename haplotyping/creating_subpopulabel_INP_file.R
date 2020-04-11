# Objective: create subpopuation labels INP file for fastPHASE input
date ()
Sys.info ()[c("nodename", "user")]
commandArgs ()
rm (list = ls ())
R.version.string ## "R version 3.6.1 (2019-07-05)"

require(dplyr); packageDescription ("dplyr", fields = "Version") #"0.8.3"

# Working directory
setwd("~/Documents/STRs/ANALYSIS/haplotyping/AR/")

# Create the dataframe we will use to name or label each population
superpopulations = c("AFR","AFR","AFR","AFR","AFR","AFR","AFR", 
                     "AMR", "AMR","AMR","AMR",
                     "EUR","EUR","EUR","EUR","EUR", 
                     "EAS","EAS","EAS","EAS","EAS",
                     "SAS","SAS","SAS","SAS","SAS")
sub_populations = c("ESN", "YRI", "GWD", "LWK", "MSL", "ACB", "ASW", 
                    "MXL", "PUR", "PEL", "CLM",
                    "IBS", "TSI", "GBR", "CEU", "FIN",
                    "JPT", "CHS", "CHB", "CDX", "KHV",
                    "PJL", "STU", "BEB", "ITU", "GIH")

# Create integer labels for each SUPER-POPULATION
superpopulations_labels = c(1,1,1,1,1,1,1,
                            2,2,2,2,
                            3,3,3,3,3,
                            4,4,4,4,4,
                            5,5,5,5,5)
# Create integer labels for each SUB-POPULATION
sub_populations_labels = c(1,2,3,4,5,6,7,
                           21,22,23,24,
                           31,32,33,34,35,
                           41,42,43,44,45,
                           51,52,53,54,55)

popu_1kg = data.frame(cbind(superpopulations, 
                            sub_populations,
                            superpopulations_labels,
                            sub_populations_labels),
                      stringsAsFactors = F)

popu_1kg$sub_populations_labels = as.integer(popu_1kg$sub_populations_labels)
popu_1kg$superpopulations_labels = as.integer(popu_1kg$superpopulations_labels)

# Load list of platekeys within VCF fo run through fastPHASE
list_platekeys = read.table("./list_33714_platekeys.txt", stringsAsFactors = F)
list_platekeys = list_platekeys$V1
length(list_platekeys)
# 33714