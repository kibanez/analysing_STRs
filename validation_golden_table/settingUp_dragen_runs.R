# Objective: setup the input file with the CRAM files we do have frmo dragen v3.2
# All have been realigned to GRCh38, in CRAM format
# The format we want to have is
# <PLATEKEY> \t <PATH> \t <GENDER> \t <GENOME_BUILD>
date ()
Sys.info ()[c("nodename", "user")]
commandArgs ()
rm (list = ls ())
R.version.string ## "R version 3.6.1 (2019-07-05)"

# libraries
library(dplyr)

# working directory
setwd("~/Documents/STRs/VALIDATION/dragen/input")

all_data = read.csv("list_fam_platekey_path_dragenv3.2.csv",
                    sep = ",",
                    stringsAsFactors = F,
                    header = F)

clin_data = read.csv("~/Documents/STRs/clinical_data/clinical_data/rd_genomes_all_data_041219.tsv",
                     sep = "\t",
                     header = T,
                     stringsAsFactors = F)

