# Objective: analyse how many family-cases we've got inside ~11k participants
# a familiar case means when anyone within the family-history (regardless of whether it's been recruited or not) is affected
# for appeal work - july 2021
date ()
Sys.info ()[c("nodename", "user")]
commandArgs ()
rm (list = ls ())
R.version.string ## "R version 3.3.2 (2016-10-31)"

# libraries
library(dplyr); packageDescription ("dplyr", fields = "Version") #"0.8.5"
library(tidyverse); packageDescription ("tidyverse", fields = "Version") # "1.2.1
library(ggplot2); packageDescription ("ggplot2", fields = "Version") #"3.3.0"

# set the working directory
setwd("~/Documents/STRs/PAPERS/VALIDATION_PAPER/")

# load dataset
l_neuro = read.table("./list_11631_PIDs.txt", stringsAsFactors = F)
l_neuro = l_neuro$V1
length(l_neuro)
# 11631