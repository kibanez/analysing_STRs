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

# Load data from every release from V2
# V2
rd_v2 = read.csv("./RE_clinical_data_V2/sequencing_report_2020-07-07_11-23-59.tsv",
                 stringsAsFactors = F,
                 header = T,
                 sep = "\t")
dim(rd_v2)
# 32436 11

rd_v2 = rd_v2 %>% select(participant_id, platekey, type)

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
# 32038 5
rd_v2 = rd_v2 %>% select(participant_id, platekey, rare_diseases_family_id, type, re_version)

# V3
rd_v3 = read.csv("./RE_clinical_data_V3/rare_disease_analysis_2020-07-07_11-16-47.tsv",
                 stringsAsFactors = F,
                 header = T,
                 sep = "\t")
dim(rd_v3)
# 51079  12
rd_v3 = rd_v3 %>% select(participant_id, rare_diseases_family_id)

cancer_v3 = read.csv("./RE_clinical_data_V3/cancer_analysis_2020-07-07_11-16-19.tsv",
                     sep = "\t",
                     stringsAsFactors = F,
                     header = T)
dim(cancer_v3)
# 15778  10

cancer_v3 = cancer_v3 %>% select(participant_id,Type)
colnames(cancer_v3) = colnames(rd_v3)
rd_v3 = rbind(rd_v3,
              cancer_v3)

seq_report_v3 = read.csv("./RE_clinical_data_V3/sequencing_report_2020-07-07_11-18-17.tsv",
                         sep = "\t",
                         stringsAsFactors = F,
                         header = T)
dim(seq_report_v3)
# 44067  11

rd_v3 = left_join(rd_v3,
                  seq_report_v3 %>% select(participant_id, platekey, type),
                  by = "participant_id")
rd_v3$re_version = rep("RE_V3", length(rd_v3$participant_id))
rd_v3 = unique(rd_v3)
dim(rd_v3)
# 73761  5
rd_v3 = rd_v3 %>% select(participant_id, platekey, rare_diseases_family_id, type, re_version)

# V4
rd_v4 = read.csv("./RE_clinical_data_V4/rare_disease_analysis_2020-07-07_10-55-22.tsv",
                 stringsAsFactors = F,
                 header = T,
                 sep = "\t")
dim(rd_v4)
# 31972  12

rd_v4 = rd_v4 %>% select(participant_id, platekey, rare_diseases_family_id)

cancer_v4 = read.csv("./RE_clinical_data_V4/cancer_analysis_2020-07-07_10-55-01.tsv",
                     sep = "\t",
                     stringsAsFactors = F,
                     header = T)
dim(cancer_v4)
# 20144 31

cancer_v4 = cancer_v4 %>% select(participant_id,platekey,type)
colnames(cancer_v4) = colnames(rd_v4)
rd_v4 = rbind(rd_v4,
              cancer_v4)

seq_report_v4 = read.csv("./RE_clinical_data_V4/sequencing_report_2020-07-07_10-56-32.tsv",
                         sep = "\t",
                         stringsAsFactors = F,
                         header = T)
dim(seq_report_v4)
# 55716  11

rd_v4 = left_join(rd_v4,
                  seq_report_v4 %>% select(participant_id, type),
                  by = "participant_id")
rd_v4$re_version = rep("RE_V4", length(rd_v4$participant_id))
rd_v4 = unique(rd_v4)
dim(rd_v4)
# 61908  5

# V5
rd_v5 = read.csv("./RE_clinical_data_V5.1/rare_disease_analysis_2020-07-07_10-47-57.tsv",
                 stringsAsFactors = F,
                 header = T,
                 sep = "\t")
dim(rd_v5)
# 68517  17

rd_v5 = rd_v5 %>% select(participant_id, platekey, rare_diseases_family_id)

cancer_v5 = read.csv("./RE_clinical_data_V5.1/cancer_analysis_2020-07-07_10-47-35.tsv",
                     sep = "\t",
                     stringsAsFactors = F,
                     header = T)
dim(cancer_v5)
# 6281  72

cancer_v5 = cancer_v5 %>% select(participant_id, germline_sample_platekey, tumour_type)
colnames(cancer_v5) = colnames(rd_v5)

rd_v5 = rbind(rd_v5,
              cancer_v5)

seq_report_v5 = read.csv("./RE_clinical_data_V5.1/sequencing_report_2020-07-07_10-48-47.tsv",
                         sep = "\t",
                         stringsAsFactors = F,
                         header = T)
dim(seq_report_v5)
# 80941  12

rd_v5 = left_join(rd_v5,
                  seq_report_v5 %>% select(participant_id, type),
                  by = "participant_id")

rd_v5$re_version = rep("RE_V5", length(rd_v5$participant_id))
rd_v5 = unique(rd_v5)
dim(rd_v5)
# 80518  5

# V6
rd_v6 = read.csv("./RE_clinical_data_V6/rare_disease_analysis_2020-07-06_16-36-31.tsv",
                 stringsAsFactors = F,
                 header = T,
                 sep = "\t")
dim(rd_v6)
# 75708  17
rd_v6 = rd_v6 %>% select(participant_id, platekey, rare_diseases_family_id)

cancer_v6 = read.csv("./RE_clinical_data_V6/cancer_analysis_2020-07-06_16-36-07.tsv",
                     sep = "\t",
                     stringsAsFactors = F,
                     header = T)
dim(cancer_v6)
# 9766  74

cancer_v6 = cancer_v6 %>% select(participant_id, germline_sample_platekey, tumour_type)
colnames(cancer_v6) = colnames(rd_v6)
rd_v6 = rbind(rd_v6,
              cancer_v6)

seq_report_v6 = read.csv("./RE_clinical_data_V6/sequencing_report_2020-07-06_16-38-35.tsv",
                         sep = "\t",
                         stringsAsFactors = F,
                         header = T)
dim(seq_report_v6)
# 91347  9

rd_v6 = left_join(rd_v6,
                  seq_report_v6 %>% select(participant_id, type),
                  by = "participant_id")
rd_v6$re_version = rep("RE_V6", length(rd_v6$participant_id))
rd_v6 = unique(rd_v6)
dim(rd_v6)
# 92958  5

# V7
rd_v7 = read.csv("./RE_clinical_data_V7/rare_disease_analysis_2020-07-06_16-26-45.tsv",
                 stringsAsFactors = F,
                 header = T,
                 sep = "\t")
dim(rd_v7)
# 74891  16

rd_v7 = rd_v7 %>% select(participant_id, plate_key, rare_diseases_family_id)

cancer_v7 = read.csv("./RE_clinical_data_V7/cancer_analysis_2020-07-06_16-26-15.tsv",
                     sep = "\t",
                     stringsAsFactors = F,
                     header = T)
dim(cancer_v7)
# 12594  77

cancer_v7 = cancer_v7 %>% select(participant_id, germline_sample_platekey, tumour_type)
colnames(cancer_v7) = colnames(rd_v7)
rd_v7 = rbind(rd_v7,
              cancer_v7)

seq_report_v7 = read.csv("./RE_clinical_data_V7/sequencing_report_2020-07-06_16-32-39.tsv",
                         sep = "\t",
                         stringsAsFactors = F,
                         header = T)
dim(seq_report_v7)
# 104907  9

rd_v7 = left_join(rd_v7,
                  seq_report_v7 %>% select(participant_id, type),
                  by = "participant_id")

rd_v7$re_version = rep("RE_V7", length(rd_v7$participant_id))
rd_v7 = unique(rd_v7)
dim(rd_v7)
# 96704  5

# Change plate_key into platekey
colnames(rd_v7) = colnames(rd_v6)

# V8
rd_v8 = read.csv("./RE_clinical_data_V8/rare_disease_analysis_2019-12-04_15-03-51.tsv",
                 stringsAsFactors = F,
                 header = T,
                 sep = "\t")
dim(rd_v8)
# 74310  16

rd_v8 = rd_v8 %>% select(participant_id, plate_key, rare_diseases_family_id)

cancer_v8 = read.csv("./RE_clinical_data_V8/cancer_analysis_2019-12-04_15-05-20.tsv",
                     sep = "\t",
                     stringsAsFactors = F,
                     header = T)
dim(cancer_v8)
# 15838  77

cancer_v8 = cancer_v8 %>% select(participant_id, germline_sample_platekey, tumour_type)
colnames(cancer_v8) = colnames(rd_v8)
rd_v8 = rbind(rd_v8,
              cancer_v8)

seq_report_v8 = read.csv("./RE_clinical_data_V8/sequencing_report_2019-12-04_15-15-48.tsv",
                         sep = "\t",
                         stringsAsFactors = F,
                         header = T)
dim(seq_report_v8)
# 107623  9

rd_v8 = left_join(rd_v8,
                  seq_report_v8 %>% select(participant_id, type),
                  by = "participant_id")
rd_v8$re_version = rep("RE_V8", length(rd_v8$participant_id))
rd_v8 = unique(rd_v8)
dim(rd_v8)
# 101623  5

colnames(rd_v8) = colnames(rd_v7)

# V9
rd_v9 = read.csv("./RE_clinical_data_V9/rare_disease_analysis_2020-07-06_15-26-22.tsv",
                 stringsAsFactors = F,
                 header = T,
                 sep = "\t")
dim(rd_v9)
# 74349  16

rd_v9 = rd_v9 %>% select(participant_id, plate_key, rare_diseases_family_id)

cancer_v9 = read.csv("./RE_clinical_data_V9/cancer_analysis_2020-07-06_15-24-44.tsv",
                     sep = "\t",
                     stringsAsFactors = F,
                     header = T)
dim(cancer_v9)
# 16347  77

cancer_v9 = cancer_v9 %>% select(participant_id, germline_sample_platekey, tumour_type)
colnames(cancer_v9) = colnames(rd_v9)
rd_v9 = rbind(rd_v9,
              cancer_v9)

seq_report_v9 = read.csv("./RE_clinical_data_V9/sequencing_report_2020-07-06_15-34-17.tsv",
                         sep = "\t",
                         stringsAsFactors = F,
                         header = T)
dim(seq_report_v9)
# 107775  9

rd_v9 = left_join(rd_v9,
                  seq_report_v9 %>% select(participant_id, type),
                  by = "participant_id")
rd_v9$re_version = rep("RE_V9", length(rd_v9$participant_id))
rd_v9 = unique(rd_v9)
dim(rd_v9)
# 102361  5

colnames(rd_v9) = colnames(rd_v7)

# V10
rd_v10 = read.csv("./RE_clinical_data_V10/rare_disease_analysis_2020-09-08_09-28-42.tsv",
                 stringsAsFactors = F,
                 header = T,
                 sep = "\t")
dim(rd_v10)
# 74271  16

rd_v10 = rd_v10 %>% select(participant_id, plate_key, rare_diseases_family_id)

cancer_v10 = read.csv("./RE_clinical_data_V10/cancer_analysis_2020-09-08_09-28-03.tsv",
                     sep = "\t",
                     stringsAsFactors = F,
                     header = T)
dim(cancer_v10)
# 16353  77

cancer_v10 = cancer_v10 %>% select(participant_id, germline_sample_platekey, tumour_type)
colnames(cancer_v10) = colnames(rd_v10)
rd_v10 = rbind(rd_v10,
              cancer_v10)

seq_report_v10 = read.csv("./RE_clinical_data_V10/sequencing_report_2020-09-08_09-31-49.tsv",
                         sep = "\t",
                         stringsAsFactors = F,
                         header = T)
dim(seq_report_v10)
# 111313  9

rd_v10 = left_join(rd_v10,
                  seq_report_v10 %>% select(participant_id, type),
                  by = "participant_id")
rd_v10$re_version = rep("RE_V10", length(rd_v10$participant_id))
rd_v10 = unique(rd_v10)
dim(rd_v10)
# 102296  5

colnames(rd_v10) = colnames(rd_v7)

# V11
rd_v11 = read.csv("./RE_clinical_data_V11/rare_disease_analysis_2020-12-30_11-25-53.tsv",
                  stringsAsFactors = F,
                  header = T,
                  sep = "\t")
dim(rd_v11)
# 74178  16

rd_v11 = rd_v11 %>% select(participant_id, plate_key, rare_diseases_family_id)

cancer_v11 = read.csv("./RE_clinical_data_V11/cancer_analysis_2020-12-30_11-25-31.tsv",
                      sep = "\t",
                      stringsAsFactors = F,
                      header = T)
dim(cancer_v11)
# 16350  77

cancer_v11 = cancer_v11 %>% select(participant_id, germline_sample_platekey, tumour_type)
colnames(cancer_v11) = colnames(rd_v11)
rd_v11 = rbind(rd_v11,
               cancer_v11)

seq_report_v11 = read.csv("./RE_clinical_data_V11/sequencing_report_2020-12-30_12-28-39.tsv",
                          sep = "\t",
                          stringsAsFactors = F,
                          header = T)
dim(seq_report_v11)
# 114238  9

rd_v11 = left_join(rd_v11,
                   seq_report_v11 %>% select(participant_id, type),
                   by = "participant_id")
rd_v11$re_version = rep("RE_V11", length(rd_v11$participant_id))
rd_v11 = unique(rd_v11)
dim(rd_v11)
# 102178  5

colnames(rd_v11) = colnames(rd_v7)

# Include here Pilot data
pilot_clin_data = read.csv("~/Documents/STRs/clinical_data/pilot_clinical_data/pilot_cohort_clinical_data_4833_genomes_removingPanels_280919.tsv",
                           stringsAsFactors = F,
                           sep = "\t",
                           header = T)
dim(pilot_clin_data)
# 4974  10

pilot_clin_data = pilot_clin_data %>% 
  select(gelID, plateKey, gelFamilyId.x)
pilot_clin_data$re_version = rep("Pilot", length(pilot_clin_data$gelID))
pilot_clin_data$type = rep("rare disease germline", length(pilot_clin_data$gelID))
pilot_clin_data = pilot_clin_data %>% select(gelID, plateKey, gelFamilyId.x, type, re_version)
colnames(pilot_clin_data) = colnames(rd_v7)

# Merge MAIN and PILOT programmes
clin_data = rbind(rd_v2,
                  rd_v3,
                  rd_v4,
                  rd_v5,
                  rd_v6,
                  rd_v7,
                  rd_v8,
                  rd_v9,
                  rd_v10,
                  rd_v11,
                  pilot_clin_data)
dim(clin_data)
# 851319  5

# Put as `,` separated the versions as pids are the same in several RE releases
list_releases = clin_data %>% 
  group_by(platekey) %>% 
  summarise(list_re_versions = toString(unique(re_version))) %>% 
  ungroup() %>% 
  as.data.frame()
dim(list_releases)
# 99534  2

clin_data = left_join(clin_data,
                      list_releases,
                      by = "platekey")

# Remove the individual `re_version`
clin_data = clin_data[,-5]

# Write into a file
write.table(clin_data,
            "~/Documents/STRs/clinical_data/clinical_data/merged_RE_releases_and_Pilot_RD_and_Cancer_PID_FID_platekey.tsv",
            sep = "\t",
            quote = F,
            row.names = F,
            col.names = T)
