# Objective: new figure 4 that can replace old Table3
date ()
Sys.info ()[c("nodename", "user")]
commandArgs ()
rm (list = ls ())
R.version.string ## "R version 3.3.2 (2016-10-31)"

# libraries
library(dplyr); packageDescription ("dplyr", fields = "Version") #"0.8.5"
library(tidyverse); packageDescription ("tidyverse", fields = "Version") # "1.2.1
library(ggplot2); packageDescription ("ggplot2", fields = "Version") #"3.3.0"
library(reshape2)
library(tidyr)

# set the working directory
setwd("~/Documents/STRs/PAPERS/VALIDATION_PAPER/LANCET/APPEAL/")
