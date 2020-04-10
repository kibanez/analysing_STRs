# Objective: create the 4 tables for the diagnostic results section
date ()
Sys.info ()[c("nodename", "user")]
commandArgs ()
rm (list = ls ())
R.version.string ## "R version 3.3.2 (2016-10-31)"

# libraries
library(dplyr); packageDescription ("dplyr", fields = "Version") #"0.8.3"
library(ggplot2); packageDescription ("ggplot2", fields = "Version") #"3.2.1"

# set the working directory
setwd("~/Documents/STRs/PAPERS/VALIDATION_PAPER/")

# Load table with the diagnostics 
# Main table
table_diseases = read.csv("table_diseases_enriched_popu_includingSkeletalMuscleChan.tsv",
                          stringsAsFactors = F, 
                          header = T,
                          sep = "\t")
dim(table_diseases)
# 11842  20

# Pilot table
table_diseases_pilot = read.csv("table_diseases_enriched_PILOT_13diseases_enriched_popu.tsv",
                                stringsAsFactors = F,
                                header = T,
                                sep = "\t")
dim(table_diseases_pilot)
# 659  13

# Define AGE, by using YOB
table_diseases = table_diseases %>%
  group_by(participant_id) %>%
  mutate(age = 2020 - year_of_birth) %>%
  ungroup() %>%
  as.data.frame()

table_diseases_pilot = table_diseases_pilot %>%
  group_by(plateKey) %>%
  mutate(age = 2020 - yearOfBirth) %>%
  ungroup() %>%
  as.data.frame()

# Define adult or paediatric
table_diseases = mutate(table_diseases, adult.paediatric = ifelse(age < 18, "Paediatric", "Adult"))
table_diseases_pilot = mutate(table_diseases_pilot, adult.paediatric = ifelse(age < 18, "Paediatric", "Adult"))

# We need to load here the repeat-size estimations for all them - EHv255
repeats_table_main = read.csv("~/Documents/STRs/data/research/EH_2.5.5_research_August2019/EH_output_v2.5.5_research_August_2019/merged_genotypeUpdated/merged_loci_86457_research_genomes_new_loci_EHv2.5.5_summer2019_removingListVcfFiles.tsv",
                              stringsAsFactors = F,
                              header = T,
                              sep = "\t")
dim(repeats_table_main)
# 3983  11

repeats_table_pilot = read.csv("~/Documents/STRs/data/pilot/EH-offtarget-v2.5.5-Pilot_October2018/merged/merged_loci_4833_Pilot_genomes_EHv2.5.5.tsv",
                               stringsAsFactors = F, 
                               header = T,
                               sep = "\t")
dim(repeats_table_pilot)
# 912  11

# Load the pathogenic threshold for the loci
gene_pathogenic_threshold = read.csv("~/git/analysing_STRs/threshold_smallest_pathogenic_reported.txt",
                                     sep = "\t",
                                     stringsAsFactors = F)

# Let's define now the 4 subtables for the purpose of the paper

# TABLE A. ONLY INCLUDING ADULTS (I.E. >= 18 IN 2020, EXCEPT FXN WHERE WE INCLUDE CHILDREN), USING FULL-MUTATION CUTOFF THRESHOLD  
# (OR YOU CAN PRODUCE A TABLE USING THE PREMUTATION CUTOFF, BUT I SUSPOECT IT WILL BE VERY NOISY AND WOULD NOT REFELCT THE THRESHOLDS THAT PANELAPP IS CURRENTLY USING)

# select diseases we are interested for TABLE A
table_a = table_diseases %>%
  filter(normalised_specific_disease %in% c("Amyotrophic lateral sclerosis or motor neuron disease", 
                                            "Charcot-Marie-Tooth disease",
                                            "Early onset dementia", 
                                            "Early onset dystonia", 
                                            "Complex Parkinsonism", 
                                            "Hereditary ataxia", 
                                            "Hereditary spastic paraplegia",
                                            "'Early onset and familial Parkinson''s Disease'"))
dim(table_a)
# 3518  21

# Complex parkinsonism is missing here
table_a = rbind(table_a,
                table_diseases %>%
                  filter(grepl("[Cc]omplex [Pp]arkin", table_diseases$normalised_specific_disease)))
dim(table_a)
# 3659  21

# Let's define list of diseases for Table A, as we have done for the genes
l_diseases_tableA = unique(table_a$normalised_specific_disease)
length(l_diseases_tableA)
# 8

# select diseases we are interested for TABLE A - PILOT
table_a_pilot = table_diseases_pilot %>%
  filter(specificDisease %in% c("Amyotrophic lateral sclerosis/motor neuron disease",
                                "Charcot-Marie-Tooth disease",
                                "Early onset dementia (encompassing fronto-temporal dementia and prion disease)",
                                "Early onset dystonia",
                                "Complex Parkinsonism (includes pallido-pyramidal syndromes)",
                                "Hereditary ataxia",
                                "Hereditary spastic paraplegia",
                                "Early onset and familial Parkinson's Disease"))
dim(table_a_pilot)
# 418  15

# Let's define the list of genes for Table A
l_genes_tableA = c("AR_CAG", "ATN1_CAG", "ATXN1_CAG", "ATXN2_CAG", "ATXN3_CAG", "ATXN7_CAG", "CACNA1A_CAG", "C9orf72_GGGGCC", "FXN_GAA", "HTT_CAG", "TBP_CAG")


# How many PIDs in the Main?
length(unique(table_a$participant_id))
# 3507

# How many PIDs are in the Pilot?
length(unique(table_a_pilot$plateKey))
# 408

# List of platekeys
# After having selected the diseases, we need to keep only with ADULTS, except for FXN we also get children -- but I'll do this a posteriori
l_platekeys_tableA = unique(table_a$plate_key.x)
length(l_platekeys_tableA)
# 3507

# PILOT
l_platekeys_tableA_pilot = unique(table_a_pilot$plateKey)
length(l_platekeys_tableA_pilot)
# 408

# Now, we want to see how many of them have an expansion on any of the genes in `l_genes_tableA`
expanded_table_main = data.frame()
for (i in 1:length(l_genes_tableA)){
  locus_name = l_genes_tableA[i]
  patho_cutoff = gene_pathogenic_threshold %>% 
    filter(locus %in% locus_name) %>%
    select(threshold) %>%
    pull()
  
  print(locus_name)
  print(patho_cutoff)
  
  expanded_table_main = rbind(expanded_table_main,
                              repeats_table_main %>% 
                                filter(gene %in% locus_name, allele >= patho_cutoff) %>%
                                select(gene, allele, Repeat_Motif, num_samples, list_samples))
  
}
dim(expanded_table_main)
# 310  5

# Now, we want to see how many of them have an expansion on any of the genes in `l_genes_tableA` - but for Pilot data
expanded_table_pilot = data.frame()
for (i in 1:length(l_genes_tableA)){
  locus_name = l_genes_tableA[i]
  patho_cutoff = gene_pathogenic_threshold %>% 
    filter(locus %in% locus_name) %>%
    select(threshold) %>%
    pull()
  
  print(locus_name)
  print(patho_cutoff)
  
  expanded_table_pilot = rbind(expanded_table_pilot,
                              repeats_table_pilot %>% 
                                filter(gene %in% locus_name, allele >= patho_cutoff) %>%
                                select(gene, allele, Repeat_Motif, num_samples, list_samples))
  
}
dim(expanded_table_pilot)
# 48  5

# After having selected the diseases, we need to keep only with ADULTS, except for FXN we also get children -- but I'll do this a posteriori
# And also, focus only in the list of platekeys of Table A
expanded_table_main_per_locus = data.frame()
index_kutre = 1
for (i in 1:length(expanded_table_main$gene)){
  list_affected_vcf = strsplit(expanded_table_main$list_samples[i], ';')[[1]]
  for (j in 1:length(list_affected_vcf)){
    expanded_table_main_per_locus = rbind(expanded_table_main_per_locus,
                                         expanded_table_main[i,])
    expanded_table_main_per_locus$list_samples[index_kutre] = sub(".vcf", "", sub("EH_", "", list_affected_vcf[j]))
    index_kutre = index_kutre + 1
  }
}
expanded_table_main_per_locus = unique(expanded_table_main_per_locus)
dim(expanded_table_main_per_locus)
# 1571  5

# The same for PILOT
expanded_table_pilot_per_locus = data.frame()
index_kutre = 1
for (i in 1:length(expanded_table_pilot$gene)){
  list_affected_vcf = strsplit(expanded_table_pilot$list_samples[i], ';')[[1]]
  for (j in 1:length(list_affected_vcf)){
    expanded_table_pilot_per_locus = rbind(expanded_table_pilot_per_locus,
                                           expanded_table_pilot[i,])
    expanded_table_pilot_per_locus$list_samples[index_kutre] = sub(".vcf", "", sub("EH_", "", list_affected_vcf[j]))
    index_kutre = index_kutre + 1
  }
}
expanded_table_pilot_per_locus = unique(expanded_table_pilot_per_locus)
dim(expanded_table_pilot_per_locus)
# 88  5

# From the expanded table, let's see how many are in l_platekeys_tableA
expanded_table_main_in_tableA = expanded_table_main_per_locus %>%
  filter(list_samples %in% l_platekeys_tableA)
dim(expanded_table_main_in_tableA)
# 114  5

# The same por PILOT
expanded_table_pilot_in_tableA = expanded_table_pilot_per_locus %>%
  filter(list_samples %in% l_platekeys_tableA_pilot)
dim(expanded_table_pilot_in_tableA)
# 12  5

# Let' enrich expanded TABLE A repeats with clinical data from `table_a`
table_a_expanded = left_join(expanded_table_main_in_tableA,
                    table_a,
                    by = c("list_samples" = "plate_key.x"))
dim(table_a_expanded)
# 120  25

# PILOT
table_a_pilot_expanded = left_join(expanded_table_pilot_in_tableA,
                                   table_a_pilot,
                                   by = c("list_samples" = "plateKey"))
dim(table_a_pilot_expanded)
# 12  19

# Let's filter out paediatric, and keep only ADULTS from this table, with exception for FXN (we keep all)
# We also focus on our list of genes
table_a_expanded = table_a_expanded %>%
  filter(gene %in% l_genes_tableA)
dim(table_a_expanded)
# 120  25

# Focus ONLY in adults
# FXN exception
table_a_FXN = table_a_expanded %>%
  filter(gene %in% "FXN_GAA")
dim(table_a_FXN)  
# 53  25

table_a_expanded = table_a_expanded %>%
  filter(adult.paediatric %in% "Adult")
dim(table_a_expanded)
# 107  25

table_a_expanded = rbind(table_a_expanded,
                         table_a_FXN)
table_a_expanded = unique(table_a_expanded)
dim(table_a_expanded)
# 112 225

# Simplify output TableA
table_a_expanded = table_a_expanded %>%
  select(list_samples, gene, allele, Repeat_Motif, participant_id, programme, genome_build, programme_consent_status, rare_diseases_family_id, biological_relationship_to_proband, 
         affection_status, participant_phenotypic_sex, year_of_birth, normalised_specific_disease, disease_sub_group, disease_group, family_group_type, family_medical_review_qc_state_code, 
         panel_list, best_guess_predicted_ancstry, self_reported, participant_ethnic_category, age, adult.paediatric)
colnames(table_a_expanded)[1] = "platekey" 
colnames(table_a_expanded)[3] = "repeat_size" 


write.table(table_a_expanded, "subtables/TableA_main.csv", quote = F, row.names = F, col.names = T, sep = ",")

# PILOT
# Let's filter out paediatric, and keep only ADULTS from this table, with exception for FXN (we keep all)
# We also focus on our list of genes
table_a_pilot_expanded = table_a_pilot_expanded %>%
  filter(gene %in% l_genes_tableA)
dim(table_a_pilot_expanded)
# 12  19

# Focus ONLY in adults
# FXN exception
table_a_pilot_FXN = table_a_pilot_expanded %>%
  filter(gene %in% "FXN_GAA")
dim(table_a_pilot_FXN)  
# 9  19

table_a_pilot_expanded = table_a_pilot_expanded %>%
  filter(adult.paediatric %in% "Adult")
dim(table_a_pilot_expanded)
# 12  19

table_a_pilot_expanded = rbind(table_a_pilot_expanded,
                               table_a_pilot_FXN)
table_a_pilot_expanded = unique(table_a_pilot_expanded)
dim(table_a_pilot_expanded)
# 12  19

# Simplify output PILOT TableA
table_a_pilot_expanded = table_a_pilot_expanded %>%
  select(list_samples, gene, allele, Repeat_Motif, gelID, gelFamilyId.x, sex, biological_relation_to_proband, disease_status, yearOfBirth, specificDisease, ageOfOnset,
         qc_state, panel_list, bestGUESS_sub_pop, bestGUESS_super_pop, age, adult.paediatric)
colnames(table_a_pilot_expanded)[1] = "platekey" 
colnames(table_a_pilot_expanded)[3] = "repeat_size" 


write.table(table_a_pilot_expanded, "subtables/TableA_pilot.csv", quote = F, row.names = F, col.names = T, sep = ",")


# This is the raw data for Table A - Main
# Let's do numbers for each locus and disease
matrix_to_print = matrix(ncol = 11, nrow = 8)
for(i in 1:length(l_diseases_tableA)){
  for (j in 1:length(l_genes_tableA)){
    number_to_print = table_a_expanded %>% 
      filter(normalised_specific_disease %in% l_diseases_tableA[i], gene %in% l_genes_tableA[j]) %>% 
      select(participant_id) %>% unique() %>% pull() %>% length()
    
    print(l_diseases_tableA[i])
    print(l_genes_tableA[j])
    print(number_to_print)
    matrix_to_print[i,j] = number_to_print
  }
}

rownames(matrix_to_print) = l_diseases_tableA
colnames(matrix_to_print) = l_genes_tableA

write.table(matrix_to_print, "./subtables/tableA_main_for_excel.tsv", sep = "\t", row.names = T, col.names = T, quote = F)

# This is the raw data for Table A - PILOT
# Let's do numbers for each locus and disease
matrix_to_print_pilot = matrix(ncol = 11, nrow = 8)
for(i in 1:length(l_diseases_tableA)){
  for (j in 1:length(l_genes_tableA)){
      number_to_print = table_a_pilot_expanded %>% 
      filter(specificDisease %in% l_diseases_tableA[i], gene %in% l_genes_tableA[j]) %>% 
      select(gelID) %>% unique() %>% pull() %>% length()
    
    print(l_diseases_tableA[i])
    print(l_genes_tableA[j])
    print(number_to_print)
    matrix_to_print_pilot[i,j] = number_to_print
  }
}

rownames(matrix_to_print_pilot) = l_diseases_tableA
colnames(matrix_to_print_pilot) = l_genes_tableA

write.table(matrix_to_print_pilot, "./subtables/tableA_pilot_for_excel.tsv", sep = "\t", row.names = T, col.names = T, quote = F)

################################################################################################################################################################
# TABLE B

################################################################################################################################################################
# TABLE C
# patients presenting with intellectual disability and or a neuromuscular phenotype were analysed for DMPK

# MAIN
table_c = table_diseases %>%
  filter(normalised_specific_disease %in% c("Intellectual disability",
                                            "Congenital muscular dystrophy",
                                            "Congenital myopathy",
                                            "Skeletal Muscle Channelopathies",
                                            "Distal myopathies"))
dim(table_c)
# 7665  21

# PILOT
table_c_pilot = table_diseases_pilot %>%
  filter(specificDisease %in% c("Intellectual disability",
                                "Congenital muscular dystrophy",
                                "Congenital myopathy",
                                "Skeletal Muscle Channelopathies",
                                "Distal myopathies"))
dim(table_c_pilot)
# 241  15

# Let's define the list of genes for Table C
l_genes_tableC = c("DMPK_CTG")

# How many PIDs in the Main?
length(unique(table_c$participant_id))
# 7317

# How many PIDs are in the Pilot?
length(unique(table_c_pilot$plateKey))
# 240 

l_platekeys_tableC = unique(table_c$plate_key.x)
l_platekeys_tableC_pilot = unique(table_c_pilot$plateKey)

# Now, we want to see how many of them have an expansion on any of the genes in `DMPK` (cutoff >= 50)
expanded_table_main = data.frame()
for (i in 1:length(l_genes_tableC)){
  locus_name = l_genes_tableC[i]
  patho_cutoff = gene_pathogenic_threshold %>% 
    filter(locus %in% locus_name) %>%
    select(threshold) %>%
    pull()
  
  print(locus_name)
  print(patho_cutoff)
  
  expanded_table_main = rbind(expanded_table_main,
                              repeats_table_main %>% 
                                filter(gene %in% locus_name, allele >= patho_cutoff) %>%
                                select(gene, allele, Repeat_Motif, num_samples, list_samples))
  
}
dim(expanded_table_main)
# 55  5

# Now, we want to see how many of them have an expansion on any of the genes in `l_genes_tableA` - but for Pilot data
expanded_table_pilot = data.frame()
for (i in 1:length(l_genes_tableC)){
  locus_name = l_genes_tableC[i]
  patho_cutoff = gene_pathogenic_threshold %>% 
    filter(locus %in% locus_name) %>%
    select(threshold) %>%
    pull()
  
  print(locus_name)
  print(patho_cutoff)
  
  expanded_table_pilot = rbind(expanded_table_pilot,
                               repeats_table_pilot %>% 
                                 filter(gene %in% locus_name, allele >= patho_cutoff) %>%
                                 select(gene, allele, Repeat_Motif, num_samples, list_samples))
  
}
dim(expanded_table_pilot)
# 2  5

# After having selected the diseases, we need to keep only with ADULTS, except for FXN we also get children -- but I'll do this a posteriori
# And also, focus only in the list of platekeys of Table A
expanded_table_main_per_locus = data.frame()
index_kutre = 1
for (i in 1:length(expanded_table_main$gene)){
  list_affected_vcf = strsplit(expanded_table_main$list_samples[i], ';')[[1]]
  for (j in 1:length(list_affected_vcf)){
    expanded_table_main_per_locus = rbind(expanded_table_main_per_locus,
                                          expanded_table_main[i,])
    expanded_table_main_per_locus$list_samples[index_kutre] = sub(".vcf", "", sub("EH_", "", list_affected_vcf[j]))
    index_kutre = index_kutre + 1
  }
}
expanded_table_main_per_locus = unique(expanded_table_main_per_locus)
dim(expanded_table_main_per_locus)
# 88  5

# The same for PILOT
expanded_table_pilot_per_locus = data.frame()
index_kutre = 1
for (i in 1:length(expanded_table_pilot$gene)){
  list_affected_vcf = strsplit(expanded_table_pilot$list_samples[i], ';')[[1]]
  for (j in 1:length(list_affected_vcf)){
    expanded_table_pilot_per_locus = rbind(expanded_table_pilot_per_locus,
                                           expanded_table_pilot[i,])
    expanded_table_pilot_per_locus$list_samples[index_kutre] = sub(".vcf", "", sub("EH_", "", list_affected_vcf[j]))
    index_kutre = index_kutre + 1
  }
}
expanded_table_pilot_per_locus = unique(expanded_table_pilot_per_locus)
dim(expanded_table_pilot_per_locus)
# 2  5

# From the expanded table, let's see how many are in l_platekeys_tableC
expanded_table_main_in_tableC = expanded_table_main_per_locus %>%
  filter(list_samples %in% l_platekeys_tableC)
dim(expanded_table_main_in_tableC)
# 16  5

# The same por PILOT
expanded_table_pilot_in_tableC = expanded_table_pilot_per_locus %>%
  filter(list_samples %in% l_platekeys_tableC_pilot)
dim(expanded_table_pilot_in_tableC)
# 0  5

# Let' enrich expanded TABLE C repeats with clinical data from `table_a`
table_c_expanded = left_join(expanded_table_main_in_tableC,
                             table_c,
                             by = c("list_samples" = "plate_key.x"))
dim(table_c_expanded)
# 16  25

# PILOT - nothing to merge

# Simplify output TableC
table_c_expanded = table_c_expanded %>%
  select(list_samples, gene, allele, Repeat_Motif, participant_id, programme, genome_build, programme_consent_status, rare_diseases_family_id, biological_relationship_to_proband, 
         affection_status, participant_phenotypic_sex, year_of_birth, normalised_specific_disease, disease_sub_group, disease_group, family_group_type, family_medical_review_qc_state_code, 
         panel_list, best_guess_predicted_ancstry, self_reported, participant_ethnic_category, age, adult.paediatric)
colnames(table_c_expanded)[1] = "platekey" 
colnames(table_c_expanded)[3] = "repeat_size" 
write.table(table_c_expanded, "subtables/TableC_main.csv", quote = F, row.names = F, col.names = T, sep = ",")

# This is the raw data for Table C - Main
# Let's do numbers for DMPK and disease

l_diseases_tableC = unique(table_c$normalised_specific_disease)
matrix_to_print = matrix(nrow  = length(l_diseases_tableC), ncol = 1)
for(i in 1:length(l_diseases_tableC)){
  for (j in 1:length(l_genes_tableC)){
    number_to_print = table_c_expanded %>% 
      filter(normalised_specific_disease %in% l_diseases_tableC[i], gene %in% l_genes_tableC[j]) %>% 
      select(participant_id) %>% unique() %>% pull() %>% length()
    
    print(l_diseases_tableC[i])
    print(l_genes_tableC[j])
    print(number_to_print)
    matrix_to_print[i,j] = number_to_print
  }
}

rownames(matrix_to_print) = l_diseases_tableC
colnames(matrix_to_print) = l_genes_tableC

write.table(matrix_to_print, "./subtables/tableC_main_for_excel.tsv", sep = "\t", row.names = T, col.names = T, quote = F)


################################################################################################################################################################
# TABLE D