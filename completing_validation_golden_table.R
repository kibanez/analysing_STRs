# Objective: Complete our validation table together with the max CI value estimated by EH for each allele
date ()
Sys.info ()[c("nodename", "user")]
commandArgs ()
rm (list = ls ())
R.version.string ## "R version 3.3.2 (2016-10-31)"

# libraries
library(dplyr)

# set the working directory
setwd("~/Documents/STRs/VALIDATION/raw_data/")

# Load validation golden table data
val_data = read.csv("../STRVALIDATION_ALLDATA_2019-10-7_ALL_kibanez.tsv",
                    header = T,
                    sep = "\t",
                    stringsAsFactors = F)
dim(val_data)
# 639 17

# Load Pilot merged table
merged_maxCI_table_pilot = read.csv("./pilot_validation/merged/merged_validation_pilot_maxCI.tsv",
                              sep = "\t",
                              header = T,
                              stringsAsFactors = F)
dim(merged_maxCI_table_pilot)
# 292 11

merged_avg_table_pilot = read.csv("./pilot_validation/merged/merged_validation_pilot_avg.tsv",
                                  sep = "\t",
                                  header = T,
                                  stringsAsFactors = F)
dim(merged_avg_table_pilot)
# 294  11

merged_avg_table_pilot = merged_avg_table_pilot %>%
  select(gene, allele, list_samples)
colnames(merged_avg_table_pilot) = c("gene", "avg_repeat", "list_samples")
merged_avg_table_pilot = merged_avg_table_pilot %>%
  mutate(unique_id = paste(gene, avg_repeat, sep = "_"))
dim(merged_avg_table_pilot)
# 294  4

merged_maxCI_table_pilot = merged_maxCI_table_pilot %>%
  select(gene, allele, list_samples)
colnames(merged_maxCI_table_pilot) = c("gene", "maxCI_repeat", "list_samples")
merged_maxCI_table_pilot = merged_maxCI_table_pilot %>%
  mutate(unique_id = paste(gene, maxCI_repeat, sep = "_"))
dim(merged_maxCI_table_pilot)
# 292  4

merged_all_pilot = full_join(merged_avg_table_pilot,
                         merged_maxCI_table_pilot,
                         by = c("unique_id"))
dim(merged_all_pilot)
# 318  7


# Load Research merged table
merged_maxCI_table_research = read.csv("../raw_data/research_validation/merged/merged_validation_research_maxCI.tsv",
                                       sep = "\t",
                                       header = T,
                                       stringsAsFactors = F)
dim(merged_maxCI_table_research)
# 862  11

merged_maxCI_table_research = merged_maxCI_table_research %>%
  select(gene, allele, list_samples)
colnames(merged_maxCI_table_research) = c("gene", "maxCI_repeat", "list_samples")
merged_maxCI_table_research = merged_maxCI_table_research %>%
  mutate(unique_id = paste(gene, maxCI_repeat, sep = "_"))
dim(merged_maxCI_table_research)
# 862  4

merged_avg_table_research = read.csv("../raw_data/research_validation/merged/merged_validation_research_avg.tsv",
                                     sep = "\t",
                                     header = T,
                                     stringsAsFactors = F)
dim(merged_avg_table_research)
# 855  11

merged_avg_table_research = merged_avg_table_research %>%
  select(gene, allele, list_samples)
colnames(merged_avg_table_research) = c("gene", "avg_repeat", "list_samples")
merged_avg_table_research = merged_avg_table_research %>%
  mutate(unique_id = paste(gene, avg_repeat, sep = "_"))
dim(merged_avg_table_research)
# 855  4

merged_all_research = full_join(merged_avg_table_research,
                                merged_maxCI_table_research,
                                by = c("unique_id"))

merged_maxCI_table = rbind(merged_maxCI_table_pilot,
                           merged_maxCI_table_research)
dim(merged_maxCI_table)
# 1154  11


# Here should be WESSEX data (It's not research nor pilot)!!


# Let's enrich the validation golden table with the max CI value for each expansion
merged_maxCI_table = merged_maxCI_table %>%
  mutate(index = 1:n())

val_data2 = val_data
val_data2$EH_a1_maxCI = rep('NA', length(val_data$EH_a1))
val_data2$EH_a2_maxCI = rep('NA', length(val_data$EH_a2))

l_platekey = unique(val_data$LP_Number)
for (i in 1:length(l_platekey)){
  val_data_platekey = val_data %>% filter(LP_Number %in% l_platekey[i])
  l_loci = unique(val_data_platekey$loci)
  for (j in 1:length(l_loci)){
    index_table = merged_maxCI_table %>% 
      filter(grepl(l_loci[j], gene), grepl(l_platekey[i],list_samples)) %>% 
      select(index)
  }
}