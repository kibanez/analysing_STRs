# Objective: Create a heatmap with the number of confirmed repeat-expansions across loci and diseases
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

# Load raw data with the confirmed RE across diseases and genes
table3 = read.csv("./template_for_figure3.tsv",
                  stringsAsFactors = F, 
                  header = T,
                  sep = "\t")
dim(table3)
# 15  14

table3_reformat = pivot_longer(data = table3, 
                                cols = -c(1),
                                names_to = "disease", 
                                values_to = "confirmed_repeats")

table3_reformat$Spec_disease = factor(table3_reformat$Spec_disease, levels = rev(unique(table3_reformat$Spec_disease)))

ggplot(data = table3_reformat, 
       mapping = aes(x = disease,
                     y = Spec_disease,
                     fill = confirmed_repeats)) +
  scale_fill_gradient(
    name = "Cor", # changes legend title
    low = "lightgray",
    high = "red",
    na.value = "lightgray",
    limit = c(min(table3_reformat$confirmed_repeats), max(table3_reformat$confirmed_repeats)),
    space = "Lab",
    guide = "colourbar"
  ) + theme_minimal() +
  geom_tile() +
  xlab(label = "") +
  ylab(label = "")
  #scale_y_reverse()






