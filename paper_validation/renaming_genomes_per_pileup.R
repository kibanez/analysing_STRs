# Objective: rename the cases/genomes that present an expansion on GEL data with EHv2
# TP -> 60 alleles
# FP -> 10 alleles
date()
Sys.info()[c("nodename", "user")]
commandArgs()
rm(list = ls())
R.version.string # "R version 3.6.2 (2019-12-12)"

# Libraries
library(dplyr)

# Working directory
setwd("~/Documents/STRs/VALIDATION/QC_visual_inspection/pileups_june/")

# Load original or raw data
# Visual_QC_EHv2.5.5_ME.xlsx as tsv
df_pos_gel = read.csv("./list_GEL_EHv2_as_positive.txt", stringsAsFactors = F, header = F, sep = "\t")
dim(df_pos_gel)
# 66  2

# Load correspondence between GEL Platekey and Paper ID
# From Table S7
df_ids = read.csv("./correspondence_gelID_paperID.tsv", stringsAsFactors = F, sep = "\t")
dim(df_ids)
# 635  3

# Prepare at the same time input data to run python vintage script across all these
df_to_write = data.frame()

for (i in 1:length(df_pos_gel$V1)){
  prefix = paste("output_EHv255_positive_GEL/EH_", df_pos_gel$V1[i], sep = "")
  orig_vcf = paste(prefix, ".vcf", sep = "")
  orig_json = paste(prefix, ".json", sep = "")
  orig_log = paste(prefix, "_alignments_relevant_reads.log", sep = "")
  
  paperID = df_ids %>% 
    filter(gel_ID %in% df_pos_gel$V1[i], locus %in% df_pos_gel$V2[i]) %>%
    select(paper_ID) %>%
    unique() %>%
    pull()
  if (length(paperID) > 0){
    prefix_new = paste("output_EHv255_positive_GEL/", paperID, sep = "")
    new_vcf = paste(prefix_new, ".vcf", sep = "")
    new_json = paste(prefix_new, ".json", sep = "")
    new_log = paste(prefix_new, "_alignments_relevant_reads.log", sep = "")
    
    file.rename(from = orig_vcf, to = new_vcf)
    file.rename(from = orig_json, to = new_json)
    file.rename(from = orig_log, to = new_log)
    
    df_to_write = rbind(df_to_write,
                        data.frame(paperID = paperID, locus = df_pos_gel$V2[i]))
  }else{
    print("error")
  }
}


write.table(df_to_write, 
            "./input_file_to_generate_pileup_from_vintage_EHv2.csv", 
            quote = F,
            sep = ",",
            row.names = F,
            col.names = F)
