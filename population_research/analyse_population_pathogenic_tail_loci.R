# Objective: we want to analyse the populations that are enriched in the pathogenic tails for each locus
date ()
Sys.info ()[c("nodename", "user")]
commandArgs ()
rm (list = ls ())
R.version.string ## "R version 3.6.1 (2019-07-05)"

library(ggplot2); packageDescription ("ggplot2", fields = "Version") #"3.2.1"
library(grid); packageDescription ("grid", fields = "Version") #"3.6.1"
library(gridExtra); packageDescription ("gridExtra", fields = "Version") #"2.3"
library(reshape2); packageDescription ("reshape2", fields = "Version") #"1.4.3"
require(dplyr); packageDescription ("dplyr", fields = "Version") #"0.8.3"

# Set working dir
setwd("~/Documents/STRs/data/research/batch_march2020/output_EHv3.2.2/merged/")

# Output directory for plots
output_folder = "./population_pathogenic_tail/"

# load data
merged_table = read.csv("./merged_92663_genomes_EHv3.2.2.tsv",
                        sep = "\t",
                        stringsAsFactors = F, 
                        header = T)
dim(merged_table)
# 8560  12

# Load popu table we have so far
popu_table = read.csv("~/Documents/STRs/ANALYSIS/population_research/MAIN_ANCESTRY/GEL_60k_germline_dataset_fine_grained_population_assignment20200224.csv",
                      stringsAsFactors = F, 
                      sep = ",",
                      header = T)
dim(popu_table)
# 59464  36

# Load clin data, `participant_ethnic_category`
clin_data = read.csv("~/Documents/STRs/clinical_data/clinical_data/rd_genomes_all_data_300320.tsv",
                     sep = "\t",
                     stringsAsFactors = F,
                     header = T)
dim(clin_data)
# 1124633  31


# Merged GRCh37 and GRCh38 tables, recoding chr names
merged_table$chr = recode(merged_table$chr,
                "1" = "chr1",
                "2" = "chr2",
                "3" = "chr3",
                "4" = "chr4",
                "5" = "chr5",
                "6" = "chr6",
                "7" = "chr7",
                "8" = "chr8",
                "9" = "chr9",
                "10" = "chr10",
                "11" = "chr11",
                "12" = "chr12",
                "13" = "chr13",
                "14" = "chr14",
                "15" = "chr15",
                "16" = "chr16",
                "17" = "chr17",
                "18" = "chr18",
                "19" = "chr19",
                "20" = "chr20",
                "21" = "chr21",
                "22" = "chr22",
                "X" = "chrX")
merged_table = merged_table %>%
  group_by(chr, gene, allele) %>%
  mutate(total_num_samples = sum(num_samples)) %>%
  ungroup() %>%
  as.data.frame() 


# HTT
merged_table_htt = merged_table %>%
  filter(gene %in% "HTT", allele >= 40)
list_vcf_patho_htt = c()
for (i in 1:length(merged_table_htt$list_samples)){
  list_vcf_patho_htt = c(list_vcf_patho_htt,
                         strsplit(merged_table_htt$list_samples[i], ';')[[1]][1])
  
}

list_vcf_patho_htt = gsub('.vcf', '', list_vcf_patho_htt)
list_vcf_patho_htt = gsub('^EH_', '', list_vcf_patho_htt)
length(list_vcf_patho_htt)
# 19

# Enrich platekeys now with ancestry info
patho_popu = popu_table %>%
  filter(ID %in% list_vcf_patho_htt) %>%
  select(ID, best_guess_predicted_ancstry, self_reported)
dim(patho_popu)
# 5  3

patho_popu2 = clin_data %>%
  filter(plate_key %in% list_vcf_patho_htt) %>%
  select(plate_key, participant_ethnic_category) 
patho_popu2 = unique(patho_popu2)
dim(patho_popu2)
# 10 2

patho_merged = left_join(patho_popu2,
                         patho_popu,
                         by = c("plate_key" = "ID"))
write.table(patho_merged, 
            "./population_pathogenic_tail/HTT_pathogenic_tail.tsv", 
            sep = "\t",
            quote = F,
            row.names = F,
            col.names = T)

# ATN1
merged_table_atn1 = merged_table %>%
  filter(gene %in% "ATN1", allele >= 48)
list_vcf_patho_atn1 = c()
for (i in 1:length(merged_table_atn1$list_samples)){
  list_vcf_patho_atn1 = c(list_vcf_patho_atn1,
                         strsplit(merged_table_atn1$list_samples[i], ';')[[1]][1])
  
}

list_vcf_patho_atn1 = gsub('.vcf', '', list_vcf_patho_atn1)
list_vcf_patho_atn1 = gsub('^EH_', '', list_vcf_patho_atn1)
length(list_vcf_patho_atn1)
# 7

# Enrich platekeys now with ancestry info
patho_popu = popu_table %>%
  filter(ID %in% list_vcf_patho_atn1) %>%
  select(ID, best_guess_predicted_ancstry, self_reported)
dim(patho_popu)
# 5  3

patho_popu2 = clin_data %>%
  filter(plate_key %in% list_vcf_patho_atn1) %>%
  select(plate_key, participant_ethnic_category) 
patho_popu2 = unique(patho_popu2)
dim(patho_popu2)
# 4 2

patho_merged = left_join(patho_popu,
                         patho_popu2,
                         by = c("ID" = "plate_key"))
write.table(patho_merged, 
            "./population_pathogenic_tail/ATN1_pathogenic_tail.tsv", 
            sep = "\t",
            quote = F,
            row.names = F,
            col.names = T)

# ATXN10
merged_table_atxn10 = merged_table %>%
  filter(gene %in% "ATXN10", allele >= 32)
list_vcf_patho_atxn10 = c()
for (i in 1:length(merged_table_atxn10$list_samples)){
  list_vcf_patho_atxn10 = c(list_vcf_patho_atxn10,
                          strsplit(merged_table_atxn10$list_samples[i], ';')[[1]][1])
  
}

list_vcf_patho_atxn10 = gsub('.vcf', '', list_vcf_patho_atxn10)
list_vcf_patho_atxn10 = gsub('^EH_', '', list_vcf_patho_atxn10)
length(list_vcf_patho_atxn10)
# 10

# Enrich platekeys now with ancestry info
patho_popu = popu_table %>%
  filter(ID %in% list_vcf_patho_atxn10) %>%
  select(ID, best_guess_predicted_ancstry, self_reported)
dim(patho_popu)
# 6  3

patho_popu2 = clin_data %>%
  filter(plate_key %in% list_vcf_patho_atxn10) %>%
  select(plate_key, participant_ethnic_category) 
patho_popu2 = unique(patho_popu2)
dim(patho_popu2)
# 8 2

patho_merged = left_join(patho_popu2,
                         patho_popu,
                         by = c("plate_key" = "ID"))
write.table(patho_merged, 
            "./population_pathogenic_tail/ATXN10_pathogenic_tail.tsv", 
            sep = "\t",
            quote = F,
            row.names = F,
            col.names = T)


# JPH3
merged_table_jph3 = merged_table %>%
  filter(gene %in% "JPH3", allele >= 39)
list_vcf_patho_jph3 = c()
for (i in 1:length(merged_table_jph3$list_samples)){
  list_vcf_patho_jph3 = c(list_vcf_patho_jph3,
                            strsplit(merged_table_jph3$list_samples[i], ';')[[1]][1])
  
}

list_vcf_patho_jph3 = gsub('.vcf', '', list_vcf_patho_jph3)
list_vcf_patho_jph3 = gsub('^EH_', '', list_vcf_patho_jph3)
length(list_vcf_patho_jph3)
# 1

# Enrich platekeys now with ancestry info
patho_popu = popu_table %>%
  filter(ID %in% list_vcf_patho_jph3) %>%
  select(ID, best_guess_predicted_ancstry, self_reported)
dim(patho_popu)
# 0 3

patho_popu2 = clin_data %>%
  filter(plate_key %in% list_vcf_patho_jph3) %>%
  select(plate_key, participant_ethnic_category) 
patho_popu2 = unique(patho_popu2)
dim(patho_popu2)
# 1 2

write.table(patho_popu2, 
            "./population_pathogenic_tail/JPH3_pathogenic_tail.tsv", 
            sep = "\t",
            quote = F,
            row.names = F,
            col.names = T)

# NOP56
merged_table_nop56 = merged_table %>%
  filter(gene %in% "NOP56", allele >= 15)
list_vcf_patho_nop56 = c()
for (i in 1:length(merged_table_nop56$list_samples)){
  list_vcf_patho_nop56 = c(list_vcf_patho_nop56,
                          strsplit(merged_table_nop56$list_samples[i], ';')[[1]][1])
  
}

list_vcf_patho_nop56 = gsub('.vcf', '', list_vcf_patho_nop56)
list_vcf_patho_nop56 = gsub('^EH_', '', list_vcf_patho_nop56)
length(list_vcf_patho_nop56)
# 10

# Enrich platekeys now with ancestry info
patho_popu = popu_table %>%
  filter(ID %in% list_vcf_patho_nop56) %>%
  select(ID, best_guess_predicted_ancstry, self_reported)
dim(patho_popu)
# 7 3

patho_popu2 = clin_data %>%
  filter(plate_key %in% list_vcf_patho_nop56) %>%
  select(plate_key, participant_ethnic_category) 
patho_popu2 = unique(patho_popu2)
dim(patho_popu2)
# 10 2

patho_merged = left_join(patho_popu2,
                         patho_popu,
                         by = c("plate_key" = "ID"))

write.table(patho_merged, 
            "./population_pathogenic_tail/NOP56_pathogenic_tail.tsv", 
            sep = "\t",
            quote = F,
            row.names = F,
            col.names = T)

