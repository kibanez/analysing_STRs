# Objective: creating tables for each locus, in order to play with cases/controls and see possible existing signals
# For EHv3.1.2
date ()
Sys.info ()[c("nodename", "user")]
commandArgs ()
rm (list = ls ())
R.version.string ## "R version 3.6.1 (2019-07-05)"

# libraries
library(dplyr)

# Set working directory
setwd("~/Documents/STRs/ANALYSIS/cases_controls/EHv3.1.2/")

# Load data
merged_data = read.csv("/Users/kibanez/Documents/STRs/data/research/EH_3.1.2_research_October2019/merged_genotypeUpdated/merged_loci_82565_research_genomes_EHv3.1.2_october2019.tsv",
                       sep = '\t',
                       header = T,
                       stringsAsFactors = F)
dim(merged_data)
# 4495  12

# Data from RE rather than from Catalog (this clinical data has been retrieved from RE on Sept 2019)
#clin_data = read.table("~/Documents/STRs/clinical_data/clinical_data/rd_genomes_all_data_250919.tsv",
clin_data = read.table("~/Documents/STRs/clinical_data/clinical_data/rd_genomes_all_data_041219.tsv",
                       sep = "\t",
                       stringsAsFactors = FALSE, 
                       header = TRUE)
dim(clin_data)  
# 1124633  28

# Let´s put all panel names into 1 single string splitted by ','
list_panels = clin_data %>% group_by(participant_id) %>% summarise(panel_list = toString(panel_name)) %>% ungroup() %>% as.data.frame()
dim(list_panels)
# 87395  2

# Let´s put all HPO terms into 1 single string splitted by ','
list_hpos = clin_data %>% group_by(participant_id) %>% summarise(hpo_list = toString(hpo_term)) %>% ungroup() %>% as.data.frame()
dim(list_hpos)
# 87395  2

# Remove the panels and hpo columns, and include the list of panels and hpo respectively
clin_data = clin_data %>% 
  select(participant_id, platekey, rare_diseases_family_id, biological_relationship_to_proband, normalised_specific_disease, specific_disease, genome_build, disease_group, disease_sub_group, year_of_birth, participant_phenotypic_sex, programme, family_group_type, affection_status)
dim(clin_data)
# 1124633  14

clin_data = left_join(clin_data,
                      list_panels,
                      by = "participant_id")
dim(clin_data)
# 1124633  15

clin_data = left_join(clin_data,
                      list_hpos,
                      by = "participant_id")
dim(clin_data)
#  1124633   16

# Let's include the population information
popu_table = read.csv("~/Documents/STRs/ANALYSIS/population_research/population_and_super-population_definitions_across_59352_WGS_REv9_271119.tsv",
                      header = T,
                      sep = "\t",
                      stringsAsFactors = F)
dim(popu_table)
# 59356  21

clin_data = left_join(clin_data,
                      popu_table %>% select(participant_id, population),
                      by = "participant_id")
dim(clin_data)
# 1124753  17

# As output
# We want to the following output - NOTE each row is an allele (!!!)
# Family ID
# Patient ID
# LP number
# Participant type (proband / mother / brother etc)
# HTT allele
# Disease Category
# Disease Subgroup
# Panels applied
# Sex
# Age at recruitment
# Age of onset
# no of family participants
# ethnic 

l_genes = sort(unique(merged_data$gene))

for (i in 1:length(l_genes)){
  locus_data = merged_data %>% filter(gene %in% l_genes[i])
  locus_data_new = data.frame()
  
  if (dim(locus_data)[1] >0){
    # For each row, we need to split/separate in many rows as alleles
    for (j in 1:length(locus_data$chr)){
      # number of samples in which a STR-expansion has been detected
      number_samp = strsplit(locus_data$list_samples[j],";")[[1]]
      
      # Clean the name of the VCF files -- removing the full path from them, and keeping only the VCF name
      number_samp = sub("^EH_", "", number_samp)
      number_samp = sub(".vcf", "", number_samp)
      
      l_which_homo = which(grepl("x2",number_samp))
      
      for (k in 1:length(number_samp)){
        # we write for each LP/sample/participant the same line/row, but also enriching with clinical data
        new_line = locus_data[j,c(1:10)]
        new_line$list_vcf_affected = number_samp[k]
        number_samp[k] = gsub("_x2", "", number_samp[k])
        
        to_include = clin_data %>% 
          filter(platekey %in% number_samp[k]) %>% 
          select(participant_id, platekey, rare_diseases_family_id, specific_disease, disease_group, disease_sub_group, year_of_birth, participant_phenotypic_sex, biological_relationship_to_proband, affection_status, family_group_type, hpo_list, panel_list, programme, genome_build, population) %>%
          unique()
        
        if (dim(to_include)[1] <= 0){
          to_include = rep('.', dim(to_include)[2])
          to_include = as.data.frame(t(as.data.frame(to_include)))
          colnames(to_include) = c("participant_id", "platekey", "rare_diseases_family_id", "specific_disease", "disease_group", "disease_sub_group", "year_of_birth", "participant_phenotypic_sex", "biological_relationship_to_proband", "affection_status", "family_group_type", "hpo_list", "panel_list", "programme", "genome_build", "population")
        }
        new_line = cbind(new_line, to_include)
        locus_data_new = rbind(locus_data_new, new_line)
        
        # If it's an alt homo (1/1) we need to repeat the line
        if (k %in% l_which_homo){
          locus_data_new = rbind(locus_data_new, new_line)
        }
      }# k
    }#j
  }# dim(locus_data) > 0
  
  # Write all `locus_data_new` output into a file
  output_file = paste("table_STR_repeat_size_each_row_allele_EHv3.1.2", l_genes[i], sep = "_")
  output_file = paste(output_file, ".tsv" , sep = ".")
  write.table(locus_data_new, output_file, sep = "\t", quote = F, row.names = F, col.names = T)
  
  # Select interested columns
  locus_data_new = locus_data_new %>%
    select(rare_diseases_family_id, participant_id, platekey, gene, Repeat_Motif, allele, specific_disease, disease_group, disease_sub_group, year_of_birth, participant_phenotypic_sex, biological_relationship_to_proband, affection_status, family_group_type, hpo_list, panel_list, programme, genome_build, population)
  
  # Adapt column names (for better understanding)
  colnames(locus_data_new)[6] = "repeat_size"
  output_file = paste("table_STR_repeat_size_each_row_allele_EHv3.1.2", l_genes[i], sep = "_")
  output_file = paste(output_file, "_simplified.tsv" , sep = "")
  write.table(locus_data_new, output_file, sep = "\t", quote = F, row.names = F, col.names = T)
  
}

