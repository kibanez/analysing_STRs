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

# Load data from every release
# V1
rd_v1 = read.csv("RE_clinical_data_V1/registration_2020-07-13_13-14-36.tsv",
                 stringsAsFactors = F, 
                 header = T,
                 sep = "\t")
dim(rd_v1)
# 18446  3

# V2
rd_v2 = read.csv("./RE_clinical_data_V2/sequencing_report_2020-07-07_11-23-59.tsv",
                 stringsAsFactors = F,
                 header = T,
                 sep = "\t")
dim(rd_v2)
# 32436 11

rd_v2 = rd_v2 %>% select(participant_id, platekey)

parti_v2 = read.csv("./RE_clinical_data_V2/participant_2020-07-07_11-23-34.tsv",
                    stringsAsFactors = F,
                    header = T,
                    sep = "\t")
rd_v2 = left_join(rd_v2,
                  parti_v2 %>% select(participant_id, rare_diseases_family_id),
                  by = "participant_id")

rd_v2$re_version = rep("RE_V2", length(rd_v2$participant_id))
rd_v2 = unique(rd_v2)
dim(rd_v2)
# 32014 4

# V3
rd_v3 = read.csv("./RE_clinical_data_V3/rare_disease_analysis_2020-07-07_11-16-47.tsv",
                 stringsAsFactors = F,
                 header = T,
                 sep = "\t")
dim(rd_v3)
# 51079  12

rd_v3 = rd_v3 %>% select(participant_id, rare_diseases_family_id)

platekey_v3 = read.csv("./RE_clinical_data_V3/genome_file_paths_and_types_2020-07-07_11-17-23.tsv",
                       stringsAsFactors = F,
                       header = T,
                       sep = "\t")

rd_v3 = left_join(rd_v3,
                  platekey_v3 %>% select(participant_id, platekey),
                  by = "participant_id")
rd_v3 = rd_v3 %>% select(participant_id, platekey, rare_diseases_family_id)
rd_v3$re_version = rep("RE_V3", length(rd_v3$participant_id))
rd_v3 = unique(rd_v3)
dim(rd_v3)
# 50497  4

# V4
rd_v4 = read.csv("./RE_clinical_data_V4/rare_disease_analysis_2020-07-07_10-55-22.tsv",
                 stringsAsFactors = F,
                 header = T,
                 sep = "\t")
dim(rd_v4)
# 31972  12

rd_v4 = rd_v4 %>% select(participant_id, platekey, rare_diseases_family_id)
rd_v4$re_version = rep("RE_V4", length(rd_v4$participant_id))
rd_v4 = unique(rd_v4)
dim(rd_v4)
# 31397  4

# V5
rd_v5 = read.csv("./RE_clinical_data_V5.1/rare_disease_analysis_2020-07-07_10-47-57.tsv",
                 stringsAsFactors = F,
                 header = T,
                 sep = "\t")
dim(rd_v5)
# 68517  17

rd_v5 = rd_v5 %>% select(participant_id, platekey, rare_diseases_family_id)
rd_v5$re_version = rep("RE_V5", length(rd_v5$participant_id))
rd_v5 = unique(rd_v5)
dim(rd_v5)
# 68003  4