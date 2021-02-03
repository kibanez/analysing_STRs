# Objective: have a unique file which merges all RE versions (Vxx) with the most important info for RD
# PID, FAMID, Platekey, RE Version
# So, whenever we want to know when the sample was removed or is included, we can easily know this
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
setwd("~/Documents/STRs/clinical_data/clinical_data/raw/")

# 