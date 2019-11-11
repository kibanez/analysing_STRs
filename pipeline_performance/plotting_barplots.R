# Objective: from EHv2.5.5 and EHv3.1.2 tables re STR calls, plot corresponding performance barplots

date ()
Sys.info ()[c("nodename", "user")]
commandArgs ()
rm (list = ls ())
R.version.string ## "R version 3.6.1 (2019-07-05)"

# libraries
library(dplyr)
library(ggplot2)
library(reshape)

setwd("/Users/kibanez/Documents/STRs/ANALYSIS/pipeline_performance/")

table_performance = read.csv("table_summary_EHv2.5.5_and_EHv3.1.2.tsv",
                             sep = "\t",
                             header = T,
                             stringsAsFactors = F)
dim(table_performance)

# melt dataframe
melt_table_performance = melt(table_performance, id=c("test"))

melt_table_performance$variable = as.character(melt_table_performance$variable)

png("TPR_FPR_comparison_EHv2_EHv3.png")
ggplot(melt_table_performance, aes(fill = test, y = value, x = variable)) +
  geom_bar(position="dodge", stat="identity") +
  xlab("") +
  ylab("") 
dev.off()