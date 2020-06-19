# Objective: creating tables for each locus, in order to play with cases/controls and see possible existing signals
# For EHv2.5.5
date ()
Sys.info ()[c("nodename", "user")]
commandArgs ()
rm (list = ls ())
R.version.string ## "R version 3.6.1 (2019-07-05)"

# libraries
library(dplyr)

# Set working directory
setwd("~/Documents/STRs/ANALYSIS/cases_controls/batch_march/EHv255")

# Load data
merged_data = read.csv("/Users/kibanez/Documents/STRs/data/research/batch_march2020/output_EHv2.5.5_vcfs/merged/merged_loci_92665_research_genomes_EHv2.5.5_batch_march2020.tsv",
                       sep = '\t',
                       header = T,
                       stringsAsFactors = F)
dim(merged_data)
# 7857  11

# Data from RE rather than from Catalog (this clinical data has been retrieved from RE on Sept 2019)
clin_data = read.table("~/Documents/STRs/clinical_data/clinical_data/rd_genomes_all_data_300320.tsv",
                       sep = "\t",
                       stringsAsFactors = FALSE, 
                       header = TRUE)
dim(clin_data)  
# 1124633  31


# Let´s put all panel names into 1 single string splitted by ','
list_panels = clin_data %>% group_by(participant_id) %>% summarise(panel_list = toString(unique(panel_name))) %>% ungroup() %>% as.data.frame()
dim(list_panels)
# 87395  2

# Let´s put all HPO terms into 1 single string splitted by ','
list_hpos = clin_data %>% group_by(participant_id) %>% summarise(hpo_list = toString(unique(hpo_term))) %>% ungroup() %>% as.data.frame()
dim(list_hpos)
# 87395  2

# Let's put all specific diseases into 1 single string splitted by ','
list_diseases = clin_data %>% group_by(participant_id) %>% summarise(diseases_list = toString(unique(normalised_specific_disease))) %>% ungroup() %>% as.data.frame()
dim(list_diseases)
# 87395  2

list_disease_group = clin_data %>% group_by(participant_id) %>% summarise(diseasegroup_list = toString(unique(disease_group))) %>% ungroup() %>% as.data.frame()
dim(list_disease_group)
# 87395  2

list_disease_subgroup = clin_data %>% group_by(participant_id) %>% summarise(diseasesubgroup_list = toString(unique(disease_sub_group))) %>% ungroup() %>% as.data.frame()
dim(list_disease_subgroup)
# 87395  2

# Remove the panels and hpo columns, and include the list of panels and hpo respectively
clin_data = clin_data %>% 
  select(participant_id, platekey, rare_diseases_family_id, biological_relationship_to_proband, genome_build, year_of_birth, participant_phenotypic_sex, programme, family_group_type, affection_status, best_guess_predicted_ancstry, self_reported)
dim(clin_data)
# 1124633  12

clin_data = left_join(clin_data,
                      list_diseases,
                      by = "participant_id")
dim(clin_data)
# 1124633  13

clin_data = left_join(clin_data,
                      list_disease_group,
                      by = "participant_id")
dim(clin_data)
# 1124633  14

clin_data = left_join(clin_data,
                      list_disease_subgroup,
                      by = "participant_id")
dim(clin_data)
# 1124633  15


clin_data = left_join(clin_data,
                      list_panels,
                      by = "participant_id")
dim(clin_data)
# 1124633  16

clin_data = left_join(clin_data,
                      list_hpos,
                      by = "participant_id")
dim(clin_data)
#  1124633   17


# Let's include the population information
popu_table = read.csv("~/Documents/STRs/ANALYSIS/population_research/MAIN_ANCESTRY/GEL_60k_germline_dataset_fine_grained_population_assignment20200224.csv",
                      header = T,
                      sep = ",",
                      stringsAsFactors = F)
dim(popu_table)
# 59464  36

clin_data = left_join(clin_data,
                      popu_table %>% select(ID, best_guess_predicted_ancstry),
                      by = c("platekey" = "ID"))
dim(clin_data)
# 1124633  18

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
        number_samp[k] = gsub("_x2", "", number_samp[k])
        new_line$list_vcf_affected = number_samp[k]
        
        to_include = clin_data %>% 
          filter(platekey %in% number_samp[k]) %>% 
          select(participant_id, platekey, rare_diseases_family_id, diseases_list, diseasegroup_list, diseasesubgroup_list, year_of_birth, participant_phenotypic_sex, biological_relationship_to_proband, affection_status, family_group_type, hpo_list, panel_list, programme, genome_build, best_guess_predicted_ancstry.x) %>%
          unique()
        
        if (dim(to_include)[1] <= 0){
          to_include = rep('.', dim(to_include)[2])
          to_include = as.data.frame(t(as.data.frame(to_include)), stringsAsFactors = F)
          colnames(to_include) = c("participant_id", "platekey", "rare_diseases_family_id", "diseases_list", "diseasegroup_list", "diseasesubgroup_list", "year_of_birth", "participant_phenotypic_sex", "biological_relationship_to_proband", "affection_status", "family_group_type", "hpo_list", "panel_list", "programme", "genome_build", "best_guess_predicted_ancstry.x")
          
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
  output_file = paste("table_STR_repeat_size_each_row_allele_EHv2.5.5", l_genes[i], sep = "_")
  output_file = paste(output_file, ".tsv" , sep = ".")
  write.table(locus_data_new, output_file, sep = "\t", quote = F, row.names = F, col.names = T)
  
  # Select interested columns
  locus_data_new = locus_data_new %>%
    select(rare_diseases_family_id, participant_id, list_vcf_affected, gene, Repeat_Motif, allele, diseases_list, diseasegroup_list, diseasesubgroup_list, year_of_birth, participant_phenotypic_sex, biological_relationship_to_proband, affection_status, family_group_type, hpo_list, panel_list, programme, genome_build, best_guess_predicted_ancstry.x)
  
  # Adapt column names (for better understanding)
  colnames(locus_data_new)[3] = "platekey"
  colnames(locus_data_new)[6] = "repeat_size"
  output_file = paste("table_STR_repeat_size_each_row_allele_EHv2.5.5", l_genes[i], sep = "_")
  output_file = paste(output_file, "_simplified.tsv" , sep = "")
  write.table(locus_data_new, output_file, sep = "\t", quote = F, row.names = F, col.names = T)
  
}




