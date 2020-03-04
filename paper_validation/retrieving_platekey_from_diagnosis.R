# Objective: retrieve the corresponding platekeys for the pids in the diagnosis table, and see the delivery version of the genomes
date ()
Sys.info ()[c("nodename", "user")]
commandArgs ()
rm (list = ls ())
R.version.string ## "R version 3.6.1 (2019-07-05)"

# libraries
library(dplyr)

# setup working directory
setwd("~/Documents/STRs/PAPERS/VALIDATION_PAPER/raw/")

# load version and genomes data
upload_report = read.csv("~/Documents/STRs/data/research/input/upload_report.020320.txt",
                         sep = "\t",
                         stringsAsFactors = F,
                         header = T)
dim(upload_report)
# 120648  10

# ACHTUNG!! this table contains everything!! there are duplicated platekeys that have been realigned to GRCh38 a posteriori --> let's take the latest one!!
upload_report = upload_report %>%
  group_by(Platekey) %>%
  mutate(latest_delivery_version = max(Delivery.Version)) %>%
  ungroup() %>%
  select(Platekey, latest_delivery_version) %>%
  as.data.frame()

upload_report = unique(upload_report)
dim(upload_report)
# 115014  2


# load RD catalog data
catalog_b38 = read.csv("~/Documents/STRs/data/research/input/batch_march2020_EHv2.5.5_and_EHv3.2.2/output_catalog_RDb38_280220.tsv",
                       sep = "\t",
                       stringsAsFactors = F,
                       header = F)
dim(catalog_b38)
# 76949  8

colnames(catalog_b38) = c("cohort_id", "platekey", "participant_id", "isProband","sex", "affection_status", "build", "programme")

# remove those cohorts having REVOKED ord DEPRECATED
catalog_b38 = catalog_b38 %>% filter(!grepl("REVOKED", cohort_id))
catalog_b38 = catalog_b38 %>% filter(!grepl("DEPRECATED", cohort_id))

# Input list of pids in the diagnosis table
list_pids = read.table("./list_pid_diagnosis.txt", stringsAsFactors = F)
list_pids = list_pids$V1
length(list_pids)
# 71

df_platekeys_pid = catalog_b38 %>% 
  filter(participant_id %in% list_pids) %>%
  select(participant_id, platekey) 
  
dim(df_platekeys_pid)
# 71  2

# Enrich them with the delivery version
df_platekeys_pid = left_join(df_platekeys_pid,
                             upload_report,
                             by = c("platekey" = "Platekey"))

write.table(df_platekeys_pid,
            "./table_pid_platekey_deliveryVersion.tsv",
            quote = F,
            sep = "\t",
            row.names = F,
            col.names = T)
