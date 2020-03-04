# Objective: retrieve the corresponding platekeys for the pids in the diagnosis table, and see the delivery version of the genomes
date ()
Sys.info ()[c("nodename", "user")]
commandArgs ()
rm (list = ls ())
R.version.string ## "R version 3.6.1 (2019-07-05)"

# libraries
library(dplyr)

# setup working directory
setwd("~/Documents/STRs/PAPERS/VALIDATION_PAPER/raw/")

upload_report = read.csv("~/Documents/STRs/data/research/input/upload_report.020320.txt",
                         sep = "\t",
                         stringsAsFactors = F,
                         header = T)
dim(upload_report)
# 120648  10

# ACHTUNG!! this table contains everything!! there are duplicated platekeys that have been realigned to GRCh38 a posteriori --> let's take the latest one!!
upload_report = upload_report %>%
  group_by(Platekey) %>%
  mutate(latest_delivery_version = max(Delivery.Version)) %>%
  ungroup() %>%
  select(Platekey, latest_delivery_version) %>%
  as.data.frame()

upload_report = unique(upload_report)
dim(upload_report)
# 115014  2
