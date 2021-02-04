# Objective: merge Pilot and Main programmes' clinical data
date()
Sys.info()[c("nodename", "user")]
commandArgs()
rm(list = ls())
R.version.string

library(dplyr)

# Set working directory
setwd("~/Documents/STRs/clinical_data/clinical_data/")

# Load latest Pilot data (frozen)


# Load latests Main clinical data release (created from our R scripts)