# Objective: rename the cases/genomes that Mike, while doing visual inspection, adapted the truth repeat-size as EH
# He changed the ‘truth’ from the PCR call
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
# List of 149 PaperIDs from Mike's work
l_ids = read.table("./list_149_GE_cases_from_Mike.txt", stringsAsFactors = F)
l_ids = l_ids$V1
length(l_ids)
# 149

# Load correspondence between GEL Platekey and Paper ID
# From Table S7
df_ids = read.csv("./correspondence_gelID_paperID.tsv", stringsAsFactors = F, sep = "\t")
dim(df_ids)
# 634  3

# Take the correspondance of gelID and locus for these 149 paper IDs
df_pos_gel = df_ids %>%
  filter(paper_ID %in% l_ids) %>%
  select(gel_ID, locus)
dim(df_pos_gel)
# 149  2

# Write down the list of gelIDs in order to bring them to local
write.table(unique(df_pos_gel$gel_ID),
            "list_unique_61_platekeys_149_GEL_platekey.txt",
            quote = F,
            row.names = F,
            col.names = F)
# Actually there are unique 61 platekeys

# Prepare at the same time input data to run python vintage script across all these
df_to_write = data.frame()

for (i in 1:length(df_pos_gel$gel_ID)){
  prefix = paste("output_list_149_Mike/EH_", df_pos_gel$gel_ID[i], sep = "")
  orig_vcf = paste(prefix, ".vcf", sep = "")
  orig_json = paste(prefix, ".json", sep = "")
  orig_log = paste(prefix, "_alignments_relevant_reads.log", sep = "")
  
  paperID = df_ids %>% 
    filter(gel_ID %in% df_pos_gel$gel_ID[i], locus %in% df_pos_gel$locus[i]) %>%
    select(paper_ID) %>%
    unique() %>%
    pull()
  if (length(paperID) > 0){
    prefix_new = paste("output_list_149_Mike/with_paperIDs/", paperID, sep = "")
    new_vcf = paste(prefix_new, ".vcf", sep = "")
    new_json = paste(prefix_new, ".json", sep = "")
    new_log = paste(prefix_new, "_alignments_relevant_reads.log", sep = "")
    
    # strategy diff, since many GE_case are same platekey but diff locus
    file.copy(from = orig_vcf, to = new_vcf)
    file.copy(from = orig_json, to = new_json)
    file.copy(from = orig_log, to = new_log)
    
    df_to_write = rbind(df_to_write,
                        data.frame(paperID = paperID, locus = df_pos_gel$locus[i]))
  }else{
    print("error")
  }
}

dim(df_to_write)
# 149  2
write.table(df_to_write, 
            "./input_file_to_generate_pileup_from_vintage_EHv2_149_Mike.csv", 
            quote = F,
            sep = ",",
            row.names = F,
            col.names = F)
