# Objective: merge into a single CSV file, the outcome of variants after being annotated with annovar and cellbase
date ()
Sys.info ()[c("nodename", "user")]
commandArgs ()
rm (list = ls ())
R.version.string ## "R version 3.6.3 (2020-02-29)"

# libraries
library(dplyr)

# defining working directory
setwd("~/Documents/STRs/ANALYSIS/FXN_monoallelic/")

# naming annovar and cellbase output files
cellbase_output = "merged_38_mono_expanded_FXN_genome_cellbase_table.csv"
annovar_output = "merged_38_mono_expanded_FXN_genome_annovar.csv.hg38_multianno.txt"

# loading annotated data 
cellbase_table = read.csv(cellbase_output,
                          header = T,
                          stringsAsFactors = F,
                          sep = ",")
dim(cellbase_table)
# 304  36

annovar_table = read.csv(annovar_output,
                         header = T,
                         stringsAsFactors = F,
                         sep = "\t")
dim(annovar_table)
# 18636  104

merged_table = left_join(cellbase_table,
                         annovar_table,
                         by = c("start" = "Start"))
dim(merged_table)
# 284 139