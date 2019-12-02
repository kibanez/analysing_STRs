# Objective: compare GEL validation golden table comparing EHv2.5.5 vs EHv3.1.2
date ()
Sys.info ()[c("nodename", "user")]
commandArgs ()
rm (list = ls ())
R.version.string ## "R version 3.6.1 (2019-07-05)"

library(ggplot2); packageDescription ("ggplot2", fields = "Version") #"2.2.1"
library(grid); packageDescription ("grid", fields = "Version") #"3.3.2"
library(gridExtra); packageDescription ("gridExtra", fields = "Version") #"2.3"
library(reshape2); packageDescription ("reshape2", fields = "Version") #"1.4.3"
require(dplyr); packageDescription ("dplyr", fields = "Version") #"0.7.4"
library(lattice); packageDescription ("lattice", fields = "Version") #"0.20-35"
library(tidyverse); packageDescription ("tidyverse", fields = "Version") #"1.2.1"

# Set working environment
setwd("/Users/kibanez/Documents/STRs/ANALYSIS/pipeline_performance/GEL_val_data_EHv255_vs_EHv312/")

# Load golden validation table - EHv2.5.5
val_data_v2 = read.csv("../EHv2_avg_VS_EHv2_maxCI/STRVALIDATION_ALLDATA_2019-10-7_ALL_kibanez_EHv255_avg_VS_EHv255_maxCI_checkFXN_withPileup_and_expValidatedData_tweaking_ATN1_updated_AR_from_NHNN.txt",
                    sep = "\t",
                    header = T,
                    stringsAsFactors = F)
dim(val_data_v2)
# 638  20

# Load golden validation table - EHv3.1.2
val_data_v3 = read.csv("../EHv3_avg_VS_EHv3_maxCI/STRVALIDATION_ALLDATA_2019-10-7_ALL_kibanez_EHv312_avg_VS_EHv312_maxCI_ClassiByAllele.tsv",
                    sep = "\t",
                    header = T,
                    stringsAsFactors = F)
dim(val_data_v3)
# 638  19

# Merge V2 and V3, after creating a new column with the EHversion
val_data_v2$EH_version = rep("EH-v2.5.5", length(val_data_v2$locus_bioinfo))
val_data_v3$EH_version = rep("EH-v3.1.2", length(val_data_v3$locus_bioinfo))

# Merge
col_intersected = intersect(colnames(val_data_v2), colnames(val_data_v3))
val_data_v2 = val_data_v2 %>% select(col_intersected)
val_data_v3 = val_data_v3 %>% select(col_intersected)

val_data = rbind(val_data_v2,
                 val_data_v3)
dim(val_data)
# 1276  16

# Filter the good ones
# 1 - Only keep with `Pileup_quality` == good or Good, `MISSING` and `blanks`
val_data_v2 = val_data %>%
  filter(Pileup_quality %in% "Good" | Pileup_quality %in% "good" | Pileup_quality %in% "" | Pileup_quality %in% "MISSING")
dim(val_data)
# 546  20

# 2 - Only keep experimental val numbers (STR_a1, STR_a2) that are integer
val_data = val_data %>%
  filter((!STR_a1 %in% "positive"))
dim(val_data)
# 538  20

val_data = val_data %>%
  filter((!STR_a2 %in% "positive"))
dim(val_data)
# 538  20

val_data = val_data %>%
  filter((!STR_a1 %in% "na"))
dim(val_data)
# 538  20

val_data = val_data %>%
  filter((!STR_a2 %in% "na"))
dim(val_data)
# 520  20

