# Objective: creating tables for each locus, in order to play with cases/controls and see possible existing signals
# For EHv3
date ()
Sys.info ()[c("nodename", "user")]
commandArgs ()
rm (list = ls ())
R.version.string ## "R version 3.6.1 (2019-07-05)"

# libraries
library(dplyr)

# Set working directory
setwd("~/Documents/STRs/ANALYSIS/cases_controls/batch_august/EHv322/")

# Load data
merged_data = read.csv("/Users/kibanez/Documents/STRs/data/research/batch_august2020/output_EHv3.2.2_vcfs/merged/merged_93446_genomes_EHv322_batch_august2020.tsv",
                       sep = '\t',
                       header = T,
                       stringsAsFactors = F)
dim(merged_data)
# 27238 12

# Data from RE rather than from Catalog (this clinical data has been retrieved from RE on Sept 2019)
# Pilot
pilot_clin_data = read.csv("~/Documents/STRs/clinical_data/pilot_clinical_data/pilot_cohort_clinical_data_4833_genomes_removingPanels_280919.tsv",
                           sep = "\t",
                           stringsAsFactors = FALSE,
                           header = TRUE)
dim(pilot_clin_data)
# 4974  10

# Let´s put all panel names into 1 single string splitted by ','
list_panels_pilot = pilot_clin_data %>% group_by(gelID) %>% summarise(panel_list = toString(unique(specificDisease))) %>% ungroup() %>% as.data.frame()
dim(list_panels_pilot)
# 4833  2

pilot_clin_data = left_join(pilot_clin_data,
                            list_panels_pilot,
                            by = "gelID")
# Remove specificDisease
pilot_clin_data = pilot_clin_data[,-8]
pilot_clin_data = unique(pilot_clin_data)
dim(pilot_clin_data)
# 4833  10

# Let's enrich with popu data
pilot_popu_table = read.csv("~/Documents/STRs/ANALYSIS/population_research/PILOT_ANCESTRY/FINE_GRAINED_RF_classifications_incl_superPOP_prediction_final20191216.csv",
                            stringsAsFactors = F,
                            header = T,
                            sep = ",")
dim(pilot_popu_table)
# 4821  44

pilot_clin_data = left_join(pilot_clin_data,
                            pilot_popu_table %>% select(ID, bestGUESS_sub_pop, bestGUESS_super_pop, self_reported),
                            by = c("plateKey"="ID"))

pilot_clin_data = unique(pilot_clin_data)
dim(pilot_clin_data)
# 4834  13

# Enrich pilot clinical data with disease group and disease subgroup
pheno_table = read.csv("~/Documents/STRs/clinical_data/pilot_clinical_data/phenotyping_v140_2019-09-13_15-26-02.tsv",
                       stringsAsFactors = F,
                       header = T,
                       sep = "\t")
pheno_table = unique(pheno_table)
dim(pheno_table)
# 106  3

pilot_clin_data = left_join(pilot_clin_data,
                            pheno_table,
                            by = c("panel_list" = "specific_disease"))
pilot_clin_data = unique(pilot_clin_data)
dim(pilot_clin_data)
# 4834  15

# Main
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

# Enrich clin_data with pilot_clin_data, keeping diff fields as `.`
colnames(pilot_clin_data) = c("participant_id", "platekey", "rare_diseases_family_id", "participant_phenotypic_sex", "biological_relationship_to_proband", "affection_status", "year_of_birth", "ageOfOnset", "qc_state", "diseases_list", "best_guess_predicted_ancstry", "bestGUESS_super_pop", "self_reported", "diseasesubgroup_list", "diseasegroup_list")

# Generate extra columns from clin data for pilot clin data
pilot_clin_data$genome_build = rep("GRCh37", length(pilot_clin_data$participant_id))
pilot_clin_data$programme = rep("RD Pilot", length(pilot_clin_data$participant_id))
pilot_clin_data$family_group_type = rep(".", length(pilot_clin_data$participant_id))
pilot_clin_data$panel_list = rep(".", length(pilot_clin_data$participant_id))
pilot_clin_data$hpo_list = rep(".", length(pilot_clin_data$participant_id))

# Select columns
pilot_clin_data = pilot_clin_data %>%
  select(participant_id, platekey, rare_diseases_family_id, biological_relationship_to_proband,
         genome_build, year_of_birth, participant_phenotypic_sex, programme,
         family_group_type, affection_status, best_guess_predicted_ancstry, self_reported,
         diseases_list, diseasegroup_list, diseasesubgroup_list, panel_list, hpo_list)

clin_data = rbind(clin_data,
                  pilot_clin_data)
dim(clin_data)
# 1129467  17

# Check genQA cases
df_genQA = read.csv("~/Documents/STRs/VALIDATION/genQA/genQA/merged_batch1_batch3_STRs.tsv",
                   stringsAsFactors = F, 
                   header = T,
                   sep = "\t")
dim(df_genQA)
# 52  17

which(df_genQA$Platekey %in% clin_data$platekey)
# there have not been included genQA cases in the batch - good!


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
          select(participant_id, platekey, rare_diseases_family_id, diseases_list, diseasegroup_list, diseasesubgroup_list, year_of_birth, participant_phenotypic_sex, biological_relationship_to_proband, affection_status, family_group_type, hpo_list, panel_list, programme, genome_build, best_guess_predicted_ancstry) %>%
          unique()
        
        if (dim(to_include)[1] <= 0){
          to_include = rep('.', dim(to_include)[2])
          to_include = as.data.frame(t(as.data.frame(to_include)), stringsAsFactors = F)
          colnames(to_include) = c("participant_id", "platekey", "rare_diseases_family_id", "diseases_list", "diseasegroup_list", "diseasesubgroup_list", "year_of_birth", "participant_phenotypic_sex", "biological_relationship_to_proband", "affection_status", "family_group_type", "hpo_list", "panel_list", "programme", "genome_build", "best_guess_predicted_ancstry")
          
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
  output_file = paste("table_STR_repeat_size_each_row_allele_EHv3.2.2", l_genes[i], sep = "_")
  output_file = paste(output_file, ".tsv" , sep = ".")
  write.table(locus_data_new, output_file, sep = "\t", quote = F, row.names = F, col.names = T)
  
  # Select interested columns
  locus_data_new = locus_data_new %>%
    select(rare_diseases_family_id, participant_id, list_vcf_affected, gene, Repeat_Motif, allele, diseases_list, diseasegroup_list, diseasesubgroup_list, year_of_birth, participant_phenotypic_sex, biological_relationship_to_proband, affection_status, family_group_type, hpo_list, panel_list, programme, genome_build, best_guess_predicted_ancstry)
  
  # Adapt column names (for better understanding)
  colnames(locus_data_new)[3] = "platekey"
  colnames(locus_data_new)[6] = "repeat_size"
  output_file = paste("table_STR_repeat_size_each_row_allele_EHv3.2.2", l_genes[i], sep = "_")
  output_file = paste(output_file, "_simplified.tsv" , sep = "")
  write.table(locus_data_new, output_file, sep = "\t", quote = F, row.names = F, col.names = T)
  
}




