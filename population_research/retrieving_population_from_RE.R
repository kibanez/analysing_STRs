# Objective: retrieve the ancestry estimation for each genome GEL research group team has estimated so far (with Plink)
date ()
Sys.info ()[c("nodename", "user")]
commandArgs ()
rm (list = ls ())
R.version.string ## "R version 3.3.2 (2016-10-31)"

# libraries
library(dplyr)

# set the working directory
setwd("~/Documents/STRs/clinical_data/clinical_data/")

popu_table = read.table("aggregate_gvcf_sample_stats_2019-10-03_22-07-31.tsv",
                        sep = "\t",
                        header = T,
                        stringsAsFactors = F)
dim(popu_table)
# 59356  51