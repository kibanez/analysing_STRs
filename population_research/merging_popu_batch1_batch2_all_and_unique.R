# Objetive: merge population batch1 and batch2, in order to have a unique space
date ()
Sys.info ()[c("nodename", "user")]
commandArgs ()
rm (list = ls ())
R.version.string ## "R version 3.6.3 (2020-02-29)"

# libraries
library(dplyr); packageDescription ("dplyr", fields = "Version") #"0.8.5"
library(ggplot2); packageDescription ("ggplot2", fields = "Version") #"3.3.0"
library(reshape); packageDescription ("reshape", fields = "Version") #"0.8.8"
library(scatterplot3d); packageDescription("scatterplot3d", fields = "Version") # 0.3-41
library(ggpubr); packageDescription("ggpubr", fields = "Version") # 0.3.0

# Set environment
setwd("/Users/kibanez/Documents/STRs/ANALYSIS/population_research/MAIN_ANCESTRY/")

# 