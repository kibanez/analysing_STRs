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
setwd("~/Documents/STRs/ANALYSIS/population_research/PAPER/expanded_genomes_main_pilot/")

# load data
merged_table = read.csv("~/Documents/STRs/data/research/batch_march2020/output_EHv3.2.2/merged/merged_92663_genomes_EHv3.2.2.tsv",
                        sep = "\t",
                        stringsAsFactors = F, 
                        header = T)
dim(merged_table)
# 8560  12

# Load MAIN popu table we have so far
popu_table = read.csv("~/Documents/STRs/ANALYSIS/population_research/MAIN_ANCESTRY/GEL_60k_germline_dataset_fine_grained_population_assignment20200224.csv",
                      stringsAsFactors = F, 
                      sep = ",",
                      header = T)
dim(popu_table)
# 59464  36

# Load PILOT popu table 
pilot_popu_table = read.csv("~/Documents/STRs/ANALYSIS/population_research/PILOT_ANCESTRY/FINE_GRAINED_RF_classifications_incl_superPOP_prediction_final20191216.csv",
                            stringsAsFactors = F,
                            sep = ",",
                            header = T)
dim(pilot_popu_table)
# 4821  44 

# Load clin data, `participant_ethnic_category`
clin_data = read.csv("~/Documents/STRs/clinical_data/clinical_data/rd_genomes_all_data_300320.tsv",
                     sep = "\t",
                     stringsAsFactors = F,
                     header = T)
dim(clin_data)
# 1124633  31

# Pilot clin data
pilot_clin_data = read.csv("~/Documents/STRs/clinical_data/pilot_clinical_data/pilot_cohort_clinical_data_4833_genomes_removingPanels_280919.tsv",
                           stringsAsFactors = F,
                           sep = "\t",
                           header = T)
dim(pilot_clin_data)
# 4974  10

# Merge popu table with family ID from clin_data
popu_table = left_join(popu_table,
                       clin_data %>% select(platekey, rare_diseases_family_id, participant_type, affection_status, normalised_specific_disease, disease_group, year_of_birth, participant_phenotypic_sex, programme, family_group_type),
                       by = c("ID" = "platekey"))
popu_table = unique(popu_table)
dim(popu_table)
# 60304  45

# For PILOT we don't have that info, but let's merge with family ID
pilot_popu_table = left_join(pilot_popu_table,
                             pilot_clin_data %>% select(plateKey, gelID, gelFamilyId.x, sex, biological_relation_to_proband, disease_status, yearOfBirth, specificDisease),
                             by = c("ID" = "plateKey"))
pilot_popu_table = unique(pilot_popu_table)
dim(pilot_popu_table)
# 4961  51

# Load unrelated list of genomes from popu - this is only main and also only the subcohort loukas' group worked on
l_unrelated = read.table("~/Documents/STRs/ANALYSIS/population_research/MAIN_ANCESTRY/60k_HWE_30k_random_unrelated_participants.txt",
                         stringsAsFactors = F)
l_unrelated = l_unrelated$V1
length(l_unrelated)
# 38344

# Let's define as unrelated: 38,344 from the main cohort + all pilot

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

popu_table = popu_table %>%
  mutate(merged_superpopu = case_when(best_guess_predicted_ancstry == "ACB" ~ "AFR",
                                      best_guess_predicted_ancstry == "ASW" ~ "AFR",
                                      best_guess_predicted_ancstry == "BEB" ~ "SAS",
                                      best_guess_predicted_ancstry == "CEU" ~ "EUR",
                                      best_guess_predicted_ancstry == "CHB" ~ "EAS",
                                      best_guess_predicted_ancstry == "CHS" ~ "EAS",
                                      best_guess_predicted_ancstry == "CLM" ~ "AMR",
                                      best_guess_predicted_ancstry == "ESN" ~ "AFR",
                                      best_guess_predicted_ancstry == "FIN" ~ "EUR",
                                      best_guess_predicted_ancstry == "GBR" ~ "EUR",
                                      best_guess_predicted_ancstry == "GIH" ~ "SAS",
                                      best_guess_predicted_ancstry == "GWD" ~ "AFR",
                                      best_guess_predicted_ancstry == "IBS" ~ "EUR",
                                      best_guess_predicted_ancstry == "ITU" ~ "SAS",
                                      best_guess_predicted_ancstry == "JPT" ~ "EAS",
                                      best_guess_predicted_ancstry == "KHV" ~ "AFR",
                                      best_guess_predicted_ancstry == "LWK" ~ "AFR",
                                      best_guess_predicted_ancstry == "MSL" ~ "AFR",
                                      best_guess_predicted_ancstry == "MXL" ~ "AMR",
                                      best_guess_predicted_ancstry == "PEL" ~ "AMR",
                                      best_guess_predicted_ancstry == "PJL" ~ "SAS",
                                      best_guess_predicted_ancstry == "PUR" ~ "AMR",
                                      best_guess_predicted_ancstry == "STU" ~ "SAS",                                      
                                      best_guess_predicted_ancstry == "TSI" ~ "EUR",
                                      best_guess_predicted_ancstry == "YRI" ~ "AFR"))


pilot_popu_table = pilot_popu_table %>%
  mutate(merged_superpopu_pilot = case_when(bestGUESS_sub_pop == "ACB" ~ "AFR",
                                      bestGUESS_sub_pop == "ASW" ~ "AFR",
                                      bestGUESS_sub_pop == "BEB" ~ "SAS",
                                      bestGUESS_sub_pop == "CEU" ~ "EUR",
                                      bestGUESS_sub_pop == "CHB" ~ "EAS",
                                      bestGUESS_sub_pop == "CHS" ~ "EAS",
                                      bestGUESS_sub_pop == "CLM" ~ "AMR",
                                      bestGUESS_sub_pop == "ESN" ~ "AFR",
                                      bestGUESS_sub_pop == "GBR" ~ "EUR",
                                      bestGUESS_sub_pop == "GIH" ~ "SAS",
                                      bestGUESS_sub_pop == "GWD" ~ "AFR",
                                      bestGUESS_sub_pop == "IBS" ~ "EUR",
                                      bestGUESS_sub_pop == "ITU" ~ "SAS",
                                      bestGUESS_sub_pop == "KHV" ~ "AFR",
                                      bestGUESS_sub_pop == "LWK" ~ "AFR",
                                      bestGUESS_sub_pop == "MSL" ~ "AFR",
                                      bestGUESS_sub_pop == "MXL" ~ "AMR",
                                      bestGUESS_sub_pop == "PJL" ~ "SAS",
                                      bestGUESS_sub_pop == "PUR" ~ "AMR",
                                      bestGUESS_sub_pop == "STU" ~ "SAS",
                                      bestGUESS_sub_pop == "TSI" ~ "EUR",
                                      bestGUESS_sub_pop == "YRI" ~ "AFR"))

# For each locus
l_genes = c("AR", "ATN1", "ATXN1", "ATXN2", "ATXN3", "ATXN7", "CACNA1A", "C9ORF72", "DMPK", "HTT", "FMR1", "FXN", "TBP")
l_premut_cutoff = c(34,34,35,31,43,17,17,30,50,35,55,44,41)

for (i in 1:length(l_genes)){
  print(l_genes[i])
  merged_table_locus = merged_table %>%
    filter(gene %in% l_genes[i], allele > l_premut_cutoff[i])
  
  list_allele_size = c()
  list_vcf_patho_locus = c()
  df_platekey_size = data.frame()
  for (j in 1:length(merged_table_locus$list_samples)){
    list_vcf_allele = strsplit(merged_table_locus$list_samples[j], ';')[[1]]
    number_vcf = length(list_vcf_allele)
    list_vcf_patho_locus = c(list_vcf_patho_locus,
                             list_vcf_allele)
    list_allele_size = rep(merged_table_locus$allele[j], number_vcf)
  
    list_vcf_allele = gsub('.vcf', '', list_vcf_allele)
    list_vcf_allele = gsub('^EH_', '', list_vcf_allele)
    list_vcf_allele = gsub('_x2', '', list_vcf_allele)
    
    # Create dataframe with platekey-repeat-size, for the expanded genomes
    df_platekey_size = rbind(df_platekey_size,
                             data.frame(platekey = list_vcf_allele,
                                  repeat_size = list_allele_size))
    df_platekey_size$platekey = as.character(df_platekey_size$platekey)
    df_platekey_size$repeat_size = as.integer(df_platekey_size$repeat_size)
  }
  
  list_vcf_patho_locus = gsub('.vcf', '', list_vcf_patho_locus)
  list_vcf_patho_locus = gsub('^EH_', '', list_vcf_patho_locus)
  list_vcf_patho_locus = gsub('_x2', '', list_vcf_patho_locus)
  
  # Enrich platekeys now with ancestry info: MAIN and PILOT
  patho_popu = popu_table %>%
    filter(ID %in% list_vcf_patho_locus) %>%
    select(ID, best_guess_predicted_ancstry, merged_superpopu, self_reported, rare_diseases_family_id, participant_type, affection_status, normalised_specific_disease, disease_group, year_of_birth, participant_phenotypic_sex, programme, family_group_type)
  print(dim(patho_popu))

  
  patho_popu2 = clin_data %>%
    filter(platekey %in% list_vcf_patho_locus) %>%
    select(platekey, participant_ethnic_category) 
  patho_popu2 = unique(patho_popu2)
  print(dim(patho_popu2))
  
  pilot_patho_popu = pilot_popu_table %>%
    filter(ID %in% list_vcf_patho_locus) %>%
    select(ID, gelID, gelFamilyId.x, sex, biological_relation_to_proband, disease_status, yearOfBirth, specificDisease, merged_superpopu_pilot, bestGUESS_sub_pop, bestGUESS_super_pop, PRED_SUM_fineGrained)
  print(dim(pilot_patho_popu))
  
  patho_merged = full_join(patho_popu,
                           patho_popu2,
                           by = c("ID" = "platekey"))
  
  patho_merged = full_join(patho_merged,
                           pilot_patho_popu,
                           by = "ID")
  
  # enrich `patho_merged` with the corresponding repeat-size for each platekey
  patho_merged = left_join(patho_merged,
                           df_platekey_size,
                           by = c("ID" = "platekey"))
  
  print(dim(patho_merged))
  
  patho_merged = unique(patho_merged)
  
  # Add locus name as column
  patho_merged$locus = rep(l_genes[i], length(patho_merged$ID))
  
  output_file_name = paste(l_genes[i], "beyond_", sep = "_")
  output_file_name = paste(output_file_name, "premutation_cutoff_", sep = "_")
  output_file_name = paste(output_file_name, as.character(l_premut_cutoff[i]), sep = "")
  output_file_name = paste(output_file_name, "EHv322_92K_population.tsv", sep = "_")
  write.table(patho_merged, 
              output_file_name, 
              sep = "\t",
              quote = F,
              row.names = F,
              col.names = T)
  
  # The same, but only keeping unrelated genomes + PILOT
  patho_merged_unrelated = patho_merged %>%
    filter(ID %in% )

  output_file_name = paste(l_genes[i], "beyond_", sep = "_")
  output_file_name = paste(output_file_name, "premutation_cutoff_", sep = "_")
  output_file_name = paste(output_file_name, as.character(l_premut_cutoff[i]), sep = "")
  output_file_name = paste(output_file_name, "EHv322_92K_population_unrelated.tsv", sep = "_")
  write.table(patho_merged_unrelated, 
              output_file_name, 
              sep = "\t",
              quote = F,
              row.names = F,
              col.names = T)
  
}


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

# Create `merged_superpopu` column 

write.table(patho_merged, 
            "./population_pathogenic_tail/NOP56_pathogenic_tail.tsv", 
            sep = "\t",
            quote = F,
            row.names = F,
            col.names = T)

# After taking the columns and merging them from Google Excel
# Set working dir
setwd("~/Documents/STRs/ANALYSIS/population_research/PAPER/expanded_genomes_main_pilot/google_excel/")

ar_table = read.csv("./batch_march_92K_EHv322_expansions_beyond_premutation - AR.csv", stringsAsFactors = F, header = T, sep = ",")
ar_table = ar_table %>% select(locus, repeat_size, ID, merged.superpopu, merged.familyID)
colnames(ar_table)[4] = "merged_superpopu"
dim(ar_table)
# 331  5

atn1_table = read.csv("./batch_march_92K_EHv322_expansions_beyond_premutation - ATN1.csv", stringsAsFactors = F, header = T, sep = ",")
atn1_table = atn1_table %>% select(locus, repeat_size, ID, merged_superpopu, merged.familyID)
dim(atn1_table)
# 30 5

atxn1_table = read.csv("./batch_march_92K_EHv322_expansions_beyond_premutation - ATXN1.csv", stringsAsFactors = F, header = T, sep = ",")
atxn1_table = atxn1_table %>% select(locus, repeat_size, ID, merged_superpopu, merged.familyID)
dim(atxn1_table)
# 5007  5

atxn2_table = read.csv("./batch_march_92K_EHv322_expansions_beyond_premutation - ATXN2.csv", stringsAsFactors = F, header = T, sep = ",")
atxn2_table = atxn2_table %>% select(locus, repeat_size, ID, merged_superpopu, merged.familyID)
dim(atxn2_table)
# 182  5

atxn3_table = read.csv("./batch_march_92K_EHv322_expansions_beyond_premutation - ATXN3.csv", stringsAsFactors = F, header = T, sep = ",")
atxn3_table = atxn3_table %>% select(locus, repeat_size, ID, merged_superpopu, merged.familyID)
dim(atxn3_table)
# 13  5

atxn7_table = read.csv("./batch_march_92K_EHv322_expansions_beyond_premutation - ATXN7.csv", stringsAsFactors = F, header = T, sep = ",")
atxn7_table = atxn7_table %>% select(locus, repeat_size, ID, merged_superpopu, merged.familyID)
dim(atxn7_table)
# 130  5

c9orf72_table = read.csv("./batch_march_92K_EHv322_expansions_beyond_premutation - C9orf72.csv", stringsAsFactors = F, header = T, sep = ",")
c9orf72_table = c9orf72_table %>% select(locus, repeat_size, ID, merged_superpopu, merged.familyID)
dim(c9orf72_table)
# 145  5

cacna1a_table = read.csv("./batch_march_92K_EHv322_expansions_beyond_premutation - CACNA1A.csv", stringsAsFactors = F, header = T, sep = ",")
cacna1a_table = cacna1a_table %>% select(locus, repeat_size, ID, merged_superpopu, merged.familyID)
dim(cacna1a_table)
# 70  5

dmpk_table = read.csv("./batch_march_92K_EHv322_expansions_beyond_premutation - DMPK.csv", stringsAsFactors = F, header = T, sep = ",")
dmpk_table = dmpk_table %>% select(locus, repeat_size, ID, merged_superpopu, merged.familyID)
dim(dmpk_table)
# 69  5

fmr1_table = read.csv("./batch_march_92K_EHv322_expansions_beyond_premutation - FMR1.csv", stringsAsFactors = F, header = T, sep = ",")
fmr1_table = fmr1_table %>% select(locus, repeat_size, ID, merged_superpopu, merged.familyID)
dim(fmr1_table)
# 736  5

fxn_table = read.csv("./batch_march_92K_EHv322_expansions_beyond_premutation - FXN.csv", stringsAsFactors = F, header = T, sep = ",")
fxn_table = fxn_table %>% select(locus, repeat_size, ID, merged_superpopu, merged.familyID)
dim(fxn_table)
# 1445  5

htt_table = read.csv("./batch_march_92K_EHv322_expansions_beyond_premutation - HTT.csv", stringsAsFactors = F, header = T, sep = ",")
htt_table = htt_table %>% select(locus, repeat_size, ID, merged_superpopu, merged.familyID)
dim(htt_table)
# 217  5

tbp_table = read.csv("./batch_march_92K_EHv322_expansions_beyond_premutation - TBP.csv", stringsAsFactors = F, header = T, sep = ",")
tbp_table = tbp_table %>% select(locus, repeat_size, ID, merged_superpopu, merged.familyID)
dim(tbp_table)
# 1276  5

merged_table = rbind(ar_table,
                     atn1_table,
                     atxn1_table,
                     atxn2_table,
                     atxn3_table,
                     atxn7_table,
                     c9orf72_table,
                     cacna1a_table,
                     dmpk_table,
                     fmr1_table,
                     fxn_table,
                     htt_table,
                     tbp_table)
dim(merged_table)
# 9651  5

# Let's focus on the ~38K unrelated genomes from Loukas' group team
l_unrelated = read.table("~/Documents/STRs/ANALYSIS/population_research/MAIN_ANCESTRY/60k_HWE_30k_random_unrelated_participants.txt", stringsAsFactors = F)
l_unrelated = l_unrelated$V1
length(l_unrelated)
# 38344

merged_table_unrelated = merged_table %>%
  filter(ID %in% l_unrelated)
dim(merged_table_unrelated)
# 4102  5

write.table(merged_table_unrelated,
            "merged_13_loci_unrelated_beyond_premutation_cutoff.tsv",
            sep = "\t",
            row.names = F,
            col.names = T)





