# Objective: we have filled in GEL genomes in golden validation table. We now want to check whether classification_a1 and classification_a2 for EHv2 and EHv3 are correct
date ()
Sys.info ()[c("nodename", "user")]
commandArgs ()
rm (list = ls ())
R.version.string ## "R version 3.6.1 (2019-07-05)"

# libraries
library(dplyr)

# set working directory
setwd("~/Documents/STRs/VALIDATION/")

# load latest validation golden table
val_data = read.csv("EHv255_EHv312_validation_cohort_GEL.tsv",
                    header = TRUE,
                    sep = "\t",
                    stringsAsFactors = FALSE)
dim(val_data)
# 638  28

# create double-check new columns in order to see whether EHv2.5.5 classifications for both alleles are OK
# create double-check new columns in order to see whether EHv3.1.2 classifications for both alleles are OK

# There are some fields not completed for `experimental_a1` and `experimental_a2`
table(val_data$experimental_a1)
# expanded   normal 
# 5       57      576 
table(val_data$experimental_a2)
# expanded   normal 
# 5        9      623 

index_to_fix = which(val_data$experimental_a1 == "")
val_data$experimental_a1[index_to_fix] = val_data$exp_PCR_a1[index_to_fix]
val_data$experimental_a2[index_to_fix] = "expanded"

val_data = val_data %>%
  group_by(LP_number, locus) %>%
  mutate(new_classification_EH2_avg1 = case_when(((EHv255_a1_avg > threshold.normal) & (experimental_a1 == "expanded")) ~ "TP",
                                             ((EHv255_a1_avg > threshold.normal) & (experimental_a1 != "expanded")) ~ "FP",
                                             ((EHv255_a1_avg < threshold.normal) & (experimental_a1 != "expanded")) ~ "TN",
                                             ((EHv255_a1_avg < threshold.normal) & (experimental_a1 == "expanded")) ~ "FN")) %>%
  as.data.frame()

val_data = val_data %>%
  group_by(LP_number, locus) %>%
  mutate(new_classification_EH2_avg2 = case_when(((EHv255_a2_avg > threshold.normal) & (experimental_a2 == "expanded")) ~ "TP",
                                             ((EHv255_a2_avg > threshold.normal) & (experimental_a2 != "expanded")) ~ "FP",
                                             ((EHv255_a2_avg < threshold.normal) & (experimental_a2 != "expanded")) ~ "TN",
                                             ((EHv255_a2_avg < threshold.normal) & (experimental_a2 == "expanded")) ~ "FN")) %>%
  as.data.frame()


val_data = val_data %>%
  group_by(LP_number, locus) %>%
  mutate(new_classification_EH3_avg1 = case_when(((EHv312_a1_avg > threshold.normal) & (experimental_a1 == "expanded")) ~ "TP",
                                                 ((EHv312_a1_avg > threshold.normal) & (experimental_a1 != "expanded")) ~ "FP",
                                                 ((EHv312_a1_avg < threshold.normal) & (experimental_a1 != "expanded")) ~ "TN",
                                                 ((EHv312_a1_avg < threshold.normal) & (experimental_a1 == "expanded")) ~ "FN")) %>%
  as.data.frame()

val_data = val_data %>%
  group_by(LP_number, locus) %>%
  mutate(new_classification_EH3_avg2 = case_when(((EHv312_a2_avg > threshold.normal) & (experimental_a2 == "expanded")) ~ "TP",
                                                 ((EHv312_a2_avg > threshold.normal) & (experimental_a2 != "expanded")) ~ "FP",
                                                 ((EHv312_a2_avg < threshold.normal) & (experimental_a2 != "expanded")) ~ "TN",
                                                 ((EHv312_a2_avg < threshold.normal) & (experimental_a2 == "expanded")) ~ "FN")) %>%
  as.data.frame()


