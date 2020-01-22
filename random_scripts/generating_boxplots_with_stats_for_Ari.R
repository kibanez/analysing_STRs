# Objective: create boxplots along with pvalue stats for case-control cases Arianna has date ()
Sys.info ()[c("nodename", "user")]
commandArgs ()
rm (list = ls ())
R.version.string ## "R version 3.6.1 (2019-07-05)"

# libraries
library(dplyr)
library(tidyverse)
library(ggplot2)
library(ggpubr)

# Set working directory
working_dir="~/Documents/STRs/ANALYSIS/cases_controls/Ari/"
setwd(working_dir)

# List all folders created by Ari so far
l_dir = list.files()

# we will go across all loci/folder
for (i in 1:length(l_dir)){
  # Let's read all case/control files
  l_files = list.files(l_dir[i])
  
  # create a new folder called `analysis`
  analysis_folder = paste(l_dir[i], "analysis", sep = "/")
  dir.create(analysis_folder)
  
  for(j in 1:length(l_files)){
    cases_motor_disorders = read.csv(paste(l_dir[i], "cases_motor_disorders_repeat_size.csv",sep = "/"), 
                                     sep = ",",
                                     stringsAsFactors = F,
                                     header = T)
    
    cases_neurodegen_motor_mitoch = read.csv(paste(l_dir[i], "cases_neurodegen_motor_mitoch_repeat_size.csv",sep = "/"), 
                                             sep = ",",
                                             stringsAsFactors = F,
                                             header = T)
    
    cases_neurodegenerative = read.csv(paste(l_dir[i], "cases_neurodegenerative_repeat_size.csv",sep = "/"), 
                                       sep = ",",
                                       stringsAsFactors = F,
                                       header = T)
    
    controls_adults_rd_cancer = read.csv(paste(l_dir[i], "controls_adults_rd_cancer_repeat_size.csv",sep = "/"), 
                                         sep = ",",
                                         stringsAsFactors = F,
                                         header = T)
    controls_adults_rd = read.csv(paste(l_dir[i], "controls_adults_rd_repeat_size.csv",sep = "/"), 
                                  sep = ",",
                                  stringsAsFactors = F,
                                  header = T)
    
    
    
    merged_table = rbind(cases_motor_disorders,
                         cases_neurodegen_motor_mitoch,
                         cases_neurodegenerative,
                         controls_adults_rd,
                         controls_adults_rd_cancer)
    
    
  }
}

