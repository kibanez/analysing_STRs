# R script to enrich with clinical data all STRs we have for all loci
# This script also generates in 1 row 1 genome or LP
# We need to split the merge TSV file by loci

# libraries
library(Rlabkey)
library(dplyr)
library(tidyverse); packageDescription ("tidyverse", fields = "Version") # "1.2.1

setwd("/Users/KristinaIbanez/Documents/STRs/GEL_STR/summer_2019/EH_2.5.5_research_August2019/EH_output_v2.5.5_research_August_2019/merged2/")

# Research 80K analysis - August 2019
merged_vcf = read.table("./merged_loci_86457_research_genomes_new_loci_EHv2.5.5_summer2019_removingListVcfFiles.tsv",
                        stringsAsFactors = F, 
                        header = T, 
                        sep = "\t")

dim(merged_vcf)                        
# 3980  11 

# Change colnames to make everything easier
colnames(merged_vcf) = c("chr", "start", "end", "repeat-size", "gene", "ref", "alt", "repeat-motif", "num_samples", "AF", "list_samples")

num_loci= unique(merged_vcf$gene)

# the Same thing, we need to stratify and distinguish between assemblies: GRCh37 and GRCh38
# This research merged TSV file is special because we do have GRCh37 and GRCh38 genomes altogether
# GRCh37 VCF files have 1,2,3,4...X,Y,MT chromosome nomenclature
# GRCh38 VCF files have chr1, chr2, chr3, ..., chrX, chrX, chrMT nomenclature
# And also, the genomic positions for the loci are different, so we cannot combine them easily

# Let's stratify the whole TSV in b37 and b38

merged_b37 = merged_vcf %>% filter(!grepl("chr", chr))
dim(merged_b37)
# 1559 11

merged_b38 = merged_vcf %>% filter(grepl("chr", chr))
dim(merged_b38)
# 2421  11


################################################################################################################################
# Data from RE rather than from Catalog (this clinical data has been retrieved from RE on Sept 2019)
clin_data = read.table("/Users/KristinaIbanez/Documents/STRs/GEL_STR/summer_2019/clinical_data/rd_genomes_all_data_250919.tsv",
                           sep = "\t",
                           stringsAsFactors = FALSE, 
                           header = TRUE)
dim(clin_data)  
# 1056568  26

# Let´s put all panel names into 1 single string splitted by ','
list_panels = clin_data %>% group_by(participant_id) %>% summarise(panel_list = toString(panel_name)) %>% ungroup() %>% as.data.frame()
dim(list_panels)
# 89163  2

# Let´s put all HPO terms into 1 single string splitted by ','
list_hpos = clin_data %>% group_by(participant_id) %>% summarise(hpo_list = toString(hpo_term)) %>% ungroup() %>% as.data.frame()
dim(list_hpos)
# 89163  2

# Remove the panels and hpo columns, and include the list of panels and hpo respectively
clin_data = clin_data %>% 
  select(participant_id, plate_key.x, rare_diseases_family_id, biological_relationship_to_proband, normalised_specific_disease, specific_disease, genome_build, disease_group, disease_sub_group, year_of_birth, participant_phenotypic_sex, programme, family_group_type, affection_status)
dim(clin_data)
# 1056568  14

clin_data = left_join(clin_data,
                      list_panels,
                      by = "participant_id")
dim(clin_data)
# 1056568  15

clin_data = left_join(clin_data,
                      list_hpos,
                      by = "participant_id")
dim(clin_data)
# 1056568  16


for (i in num_loci){
  # GRCh37 annotation
  merged_vcf = merged_b37
  
  # we go locus by locus
  df_locus = merged_vcf %>% filter(gene %in% i)
  df_locus_new = data.frame()
  
  if (dim(df_locus)[1] >0){
    
  # For each row, we need to split/separate in many rows as samples/LPs
  for (j in 1:length(df_locus$chr)){
    
    # number of samples in which a STR-expansion has been detected
    number_samp = strsplit(df_locus$list_samples[j],";")[[1]]
    
    # Clean the name of the VCF files -- removing the full path from them, and keeping only the VCF name
    number_samp = sub("^EH_", "", number_samp)
    number_samp = sub(".vcf", "", number_samp)
    
    # Include in df_locus$list_samples the number_samp
    df_locus$list_samples[j] = paste(number_samp, collapse = ';')
    
    for (k in 1:length(number_samp)){
      # we write for each LP/sample/participant the same line/row, but also enriching with clinical data
      new_line = df_locus[j,]
      new_line$list_vcf_affected = number_samp[k]
      
      # we add now clinical data for each row -- each participant/sample/genome
      
      # sometimes we do have `_x2` for the sample, if 2 alleles are having the same repeat-size (GT = 1/1)
      number_samp[k] = gsub("`_x2", "", number_samp[k])
      number_samp[k] = gsub("_x2", "", number_samp[k])
      to_include = clin_data %>% 
        filter(plate_key.x %in% number_samp[k]) %>% 
        select(participant_id, plate_key.x, rare_diseases_family_id, specific_disease, disease_group, disease_sub_group, year_of_birth, participant_phenotypic_sex, biological_relationship_to_proband, affection_status, family_group_type, hpo_list, panel_list, programme, genome_build) %>%
        unique()
    
      if (dim(to_include)[1] <= 0){
        to_include = rep('.', dim(to_include)[2])
        to_include = as.data.frame(t(as.data.frame(to_include)))
        colnames(to_include) = c("participant_id", "plate_key.x", "rare_diseases_family_id", "specific_disease", "disease_group", "disease_sub_group", "year_of_birth", "participant_phenotypic_sex", "biological_relationship_to_proband", "affection_status", "family_group_type", "hpo_list", "panel_list", "programme", "genome_build")
      }
      new_line = cbind(new_line, to_include)
      df_locus_new = rbind(df_locus_new, new_line)
    }# k
  }#j
  
  # We write for each locus the output dataframe `df_locus_new`
  output_name = paste("research_genomes_86457_GRCh37_EHv2.5.5", i, sep = "_")
  output_name = paste(output_name, "tsv", sep = ".")
  
  # we filter out the column `list_samples` from the final dataframe, not to have too much data
  df_locus_new  = df_locus_new[,-11]
  
  write.table(df_locus_new, file = output_name, quote = F, row.names = F, col.names = T, sep = "\t")
  } # if dim > 0
  ################################################################################################################################################
  ################################################################################################################################################
  
  # GRCh38 annotation
  
  merged_vcf = merged_b38
  
  # we go locus by locus
  df_locus = merged_vcf %>% filter(gene %in% i)
  df_locus_new = data.frame()
  
  if (dim(df_locus)[1] > 0){
  # For each row, we need to split/separate in many rows as samples/LPs
  for (j in 1:length(df_locus$chr)){
    
    # number of samples in which a STR-expansion has been detected
    #number_samp = strsplit(df_locus$list_vcf_affected[j],";")[[1]]
    number_samp = strsplit(df_locus$list_samples[j],";")[[1]]
    
    # Clean the name of the VCF files -- removing the full path from them, and keeping only the VCF name
    number_samp = sub("^EH_", "", number_samp)
    number_samp = sub(".vcf", "", number_samp)
    
    # Include in df_locus$list_samples the number_samp
    df_locus$list_samples[j] = paste(number_samp, collapse = ';')
    
    for (k in 1:length(number_samp)){
      # we write for each LP/sample/participant the same line/row, but also enriching with clinical data
      
      new_line = df_locus[j,]
      new_line$list_vcf_affected = number_samp[k]
      
      # we add now clinical data for each row -- each participant/sample/genome
      
      # sometimes we do have `_x2` for the sample, if 2 alleles are having the same repeat-size (GT = 1/1)
      number_samp[k] = gsub("`_x2", "", number_samp[k])
      number_samp[k] = gsub("`_x2", "", number_samp[k])
      to_include = clin_data %>% 
        filter(plate_key.x %in% number_samp[k]) %>% 
        select(participant_id, plate_key.x, rare_diseases_family_id, specific_disease, disease_group, disease_sub_group, year_of_birth, participant_phenotypic_sex, biological_relationship_to_proband, affection_status, family_group_type, hpo_list, panel_list, programme, genome_build) %>%
        unique()
      
      if (dim(to_include)[1] <= 0){
        to_include = rep('.', dim(to_include)[2])
        to_include = as.data.frame(t(as.data.frame(to_include)))
        colnames(to_include) = c("participant_id", "plate_key.x", "rare_diseases_family_id", "specific_disease", "disease_group", "disease_sub_group", "year_of_birth", "participant_phenotypic_sex", "biological_relationship_to_proband", "affection_status", "family_group_type", "hpo_list", "panel_list", "programme", "genome_build")
      }
      new_line = cbind(new_line, to_include)
      df_locus_new = rbind(df_locus_new, new_line)
      
    }# k
    
  }#j
  
  # We write for each locus the output dataframe `df_locus_new`
  output_name = paste("research_genomes_86457_GRCh38_EHv2.5.5", i, sep = "_")
  output_name = paste(output_name, "tsv", sep = ".")
  
  # we filter out the column `list_samples` from the final dataframe, not to have too much data
  df_locus_new  = df_locus_new[,-11]
  
  write.table(df_locus_new, file = output_name, quote = F, row.names = F, col.names = T, sep = "\t")
  
  }# if dim > 0
}#i







