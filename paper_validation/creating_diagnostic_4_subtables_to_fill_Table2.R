# Objective: create the 4 tables for the diagnostic results section
# Plus: enrich with number of female/male, ethnias recoded, and number of UNIQUE participants
date ()
Sys.info ()[c("nodename", "user")]
commandArgs ()
rm (list = ls ())
R.version.string ## "R version 3.3.2 (2016-10-31)"

# libraries
library(dplyr); packageDescription ("dplyr", fields = "Version") #"0.8.5"
library(ggplot2); packageDescription ("ggplot2", fields = "Version") #"3.3.0"

# set the working directory
setwd("~/Documents/STRs/PAPERS/VALIDATION_PAPER/")

# Load table with the diagnostics 
# Main table
table_diseases = read.csv("table_diseases_enriched_popu_includingSkeletalMuscleChan.tsv",
                          stringsAsFactors = F, 
                          header = T,
                          sep = "\t")
dim(table_diseases)
# 12254  19

# Pilot table
table_diseases_pilot = read.csv("table_diseases_enriched_PILOT_13diseases_enriched_popu.tsv",
                                stringsAsFactors = F,
                                header = T,
                                sep = "\t")
dim(table_diseases_pilot)
# 660  13

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

# recode repeats_table_main: FMR1 for b37 is `FMR1` and `FMR1_CGG` for b38, let's recode
repeats_table_main$gene = recode(repeats_table_main$gene,
                                 "FMR1" = "FMR1_CGG")

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


# Let's see number of participants fROM ROW DATA
# Let's simplify the table
table_a = table_a %>%
  select(participant_id, normalised_specific_disease, adult.paediatric, participant_phenotypic_sex, participant_ethnic_category, age)
table_a_pilot = table_a_pilot %>%
  select(gelID, specificDisease, adult.paediatric, sex, age)


# Let's recode the ethnicity, simplifying it
table_a$participant_ethnic_category = recode(table_a$participant_ethnic_category,
                                             "White: British"= "White",
                                             "White White: British"= "White",
                                             "White: Any other White background"="White",
                                             "Asian or Asian British: Pakistani"="Asian",
                                             "Asian or Asian British: Indian"="Asian",
                                             "Asian or Asian British: Any other Asian background"="Asian",
                                             "Black or Black British: African"="Black",
                                             "Other Ethnic Groups: Any other ethnic group"="Other", 
                                             "Mixed: Any other mixed background"="Mixed",
                                             "White: Irish"="White",
                                             "Asian or Asian British: Bangladeshi"="Asian",
                                             "Mixed: White and Asian"="Mixed",
                                             "Mixed: White and Black Caribbean"="Mixed",
                                             "Black or Black British: Caribbean"="Black",
                                             "Mixed: White and Black African"="Mixed",
                                             "Black or Black British: Any other Black background"="Black",
                                             "Other Ethnic Groups: Chinese"="Other")
# Defining NA's as `Not stated`
which_na = which(is.na(table_a$participant_ethnic_category))
table_a$participant_ethnic_category[which_na] = "Not Stated"

table_a = unique(table_a)
table_a_pilot = unique(table_a_pilot)

table_a %>% filter(normalised_specific_disease %in% "Hereditary ataxia") %>% select(normalised_specific_disease) %>% unique()
l_pid_ha_main = table_a %>% filter(normalised_specific_disease %in% "Hereditary ataxia") %>% select(participant_id) %>% unique() %>% pull() 
length(l_pid_ha_main)
# 1038
#Pilot
table_a_pilot %>% filter(specificDisease %in% "Hereditary ataxia")  %>% select(specificDisease) %>% unique()
l_pid_ha_pilot = table_a_pilot %>% filter(specificDisease %in% "Hereditary ataxia")  %>% select(gelID) %>% unique() %>% pull()
length(l_pid_ha_pilot)
# 144

# Age
# MAIN
l_age_ha_main = table_a %>% filter(participant_id %in% l_pid_ha_main) %>% select(age) %>% pull()
# 1041
l_age_ha_pilot = table_a_pilot %>% filter(gelID %in% l_pid_ha_pilot) %>% select(age) %>% pull()
# 152

l_age_merged_ha = c(l_age_ha_main,
                    l_age_ha_pilot)
length(l_age_merged_ha)

mean(l_age_merged_ha)
# 53.6
summary(l_age_merged_ha)
#  Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 5.00   41.00   57.00   53.64   70.00  101.00 

# Gender
# Main
table_a %>% filter(participant_id %in% l_pid_ha_main, participant_phenotypic_sex %in% "Female") %>% select(participant_id) %>% unique() %>% pull() %>% length()
# 513
table_a %>% filter(participant_id %in% l_pid_ha_main, participant_phenotypic_sex %in% "Male") %>% select(participant_id) %>% unique() %>% pull() %>% length()
# 525

# Pilot
table_a_pilot %>% filter(gelID %in% l_pid_ha_pilot, sex %in% "female") %>% select(gelID) %>% unique() %>% pull() %>% length()
# 72
table_a_pilot %>% filter(gelID %in% l_pid_ha_pilot, sex %in% "male") %>% select(gelID) %>% unique() %>% pull() %>% length()
# 72

# Female (main + pilot) vs total
(513 + 72) / 1182
# 0.49
# Male (main + pilot) vs total
(525 + 72) / 1182
# 0.51

# Ethnicity
# MAIN
l_eth_main = table_a %>% filter(participant_id %in% l_pid_ha_main) %>% select(participant_ethnic_category) %>% pull() #%>% table() %>% prop.table()
# PILOT - consider all them `Not stated`
l_eth_pilot = rep("Not Stated", length(table_a_pilot$gelID))

l_eth_merged = c(l_eth_main, l_eth_pilot)
prop.table(table(l_eth_merged))
#Asian       Black       Mixed  Not Stated       Other       White 
#0.052090473 0.007539411 0.008224812 0.363262509 0.006854010 0.562028787 

# Let's automate for the rest of diseases: ONlY adults (removing HA from disease list)
# ONLY adults
l_diseases_tableA = l_diseases_tableA[-2]
for(i in 1:length(l_diseases_tableA)){
  # Number of UNIQUE pids
  # MAIN
  table_a %>% filter(normalised_specific_disease %in% l_diseases_tableA[i]) %>% select(normalised_specific_disease) %>% unique() %>% print()
  l_pid_disease_main = table_a %>% filter(normalised_specific_disease %in% l_diseases_tableA[i], adult.paediatric %in% "Adult") %>% select(participant_id) %>% unique() %>% pull() 
  
  # Pilot
  table_a_pilot %>% filter(specificDisease %in% l_diseases_tableA[i])  %>% select(specificDisease) %>% unique() %>% print()
  l_pid_disease_pilot = table_a_pilot %>% filter(specificDisease %in% l_diseases_tableA[i], adult.paediatric %in% "Adult")  %>% select(gelID) %>% unique() %>% pull()
  
  # MERGE
  l_pid_merged = c(l_pid_disease_main, l_pid_disease_pilot)
  print(length(l_pid_merged))  
  
  # Age
  # MAIN
  l_age_disease_main = table_a %>% filter(participant_id %in% l_pid_disease_main) %>% select(age) %>% pull()
  l_age_disease_pilot = table_a_pilot %>% filter(gelID %in% l_pid_disease_pilot) %>% select(age) %>% pull()
  
  l_age_merged_disease = c(l_age_disease_main,
                           l_age_disease_pilot)
  print(mean(l_age_merged_disease))
  print(summary(l_age_merged_disease))
  
  # GENDER
  # Main
  female_main = table_a %>% filter(participant_id %in% l_pid_disease_main, participant_phenotypic_sex %in% "Female") %>% select(participant_id) %>% unique() %>% pull() %>% length()
  male_main = table_a %>% filter(participant_id %in% l_pid_disease_main, participant_phenotypic_sex %in% "Male") %>% select(participant_id) %>% unique() %>% pull() %>% length()
  
  # Pilot
  female_pilot = table_a_pilot %>% filter(gelID %in% l_pid_disease_pilot, sex %in% "female") %>% select(gelID) %>% unique() %>% pull() %>% length()
  male_pilot = table_a_pilot %>% filter(gelID %in% l_pid_disease_pilot, sex %in% "male") %>% select(gelID) %>% unique() %>% pull() %>% length()
  
  # Female (main + pilot) vs total
  print((female_main + female_pilot) / length(l_pid_merged))
  # Male (main + pilot) vs total
  print((male_main + male_pilot) / length(l_pid_merged))
  
  # Ethnicity
  # MAIN
  l_eth_main = table_a %>% filter(participant_id %in% l_pid_disease_main) %>% select(participant_ethnic_category) %>% pull() 
  # PILOT - consider all them `Not stated`
  l_eth_pilot = rep("Not Stated", length(table_a_pilot$gelID))
  
  l_eth_merged = c(l_eth_main, l_eth_pilot)
  print(prop.table(table(l_eth_merged)))
}


# total number for all
# HA
l_pid_ha_main = table_a %>% filter(normalised_specific_disease %in% "Hereditary ataxia") %>% select(participant_id) %>% unique() %>% pull() 
l_pid_ha_pilot = table_a_pilot %>% filter(specificDisease %in% "Hereditary ataxia")  %>% select(gelID) %>% unique() %>% pull()

# "Early onset dystonia"
l_pid_eod_main = table_a %>% filter(normalised_specific_disease %in% "Early onset dystonia", adult.paediatric %in% "Adult") %>% select(participant_id) %>% unique() %>% pull() 
l_pid_eod_pilot = table_a_pilot %>% filter(specificDisease %in% "Early onset dystonia", adult.paediatric %in% "Adult")  %>% select(gelID) %>% unique() %>% pull()

# "Hereditary spastic paraplegia"
l_pid_hsp_main = table_a %>% filter(normalised_specific_disease %in% "Hereditary spastic paraplegia", adult.paediatric %in% "Adult") %>% select(participant_id) %>% unique() %>% pull() 
l_pid_hsp_pilot = table_a_pilot %>% filter(specificDisease %in% "Hereditary spastic paraplegia", adult.paediatric %in% "Adult")  %>% select(gelID) %>% unique() %>% pull()

# "Charcot-Marie-Tooth disease" 
l_pid_cmt_main = table_a %>% filter(normalised_specific_disease %in% "Charcot-Marie-Tooth disease", adult.paediatric %in% "Adult") %>% select(participant_id) %>% unique() %>% pull() 
l_pid_cmt_pilot = table_a_pilot %>% filter(specificDisease %in% "Charcot-Marie-Tooth disease", adult.paediatric %in% "Adult")  %>% select(gelID) %>% unique() %>% pull()

# "Early onset dementia"
l_pid_eode_main = table_a %>% filter(normalised_specific_disease %in% "Early onset dementia", adult.paediatric %in% "Adult") %>% select(participant_id) %>% unique() %>% pull() 
l_pid_eode_pilot = table_a_pilot %>% filter(specificDisease %in% "Early onset dementia", adult.paediatric %in% "Adult")  %>% select(gelID) %>% unique() %>% pull()

# "Amyotrophic lateral sclerosis or motor neuron disease"
l_pid_als_main = table_a %>% filter(normalised_specific_disease %in% "Amyotrophic lateral sclerosis or motor neuron disease", adult.paediatric %in% "Adult") %>% select(participant_id) %>% unique() %>% pull() 
l_pid_als_pilot = table_a_pilot %>% filter(specificDisease %in% "Amyotrophic lateral sclerosis or motor neuron disease", adult.paediatric %in% "Adult")  %>% select(gelID) %>% unique() %>% pull()

# "'Early onset and familial Parkinson''s Disease'"
l_pid_park_main = table_a %>% filter(normalised_specific_disease %in% "'Early onset and familial Parkinson''s Disease'", adult.paediatric %in% "Adult") %>% select(participant_id) %>% unique() %>% pull() 
l_pid_park_pilot = table_a_pilot %>% filter(specificDisease %in% "'Early onset and familial Parkinson''s Disease'", adult.paediatric %in% "Adult")  %>% select(gelID) %>% unique() %>% pull()

# "Complex Parkinsonism (includes pallido-pyramidal syndromes)"
l_pid_cpark_main = table_a %>% filter(normalised_specific_disease %in% "Complex Parkinsonism (includes pallido-pyramidal syndromes)", adult.paediatric %in% "Adult") %>% select(participant_id) %>% unique() %>% pull() 
l_pid_cpark_pilot = table_a_pilot %>% filter(specificDisease %in% "Complex Parkinsonism (includes pallido-pyramidal syndromes)", adult.paediatric %in% "Adult")  %>% select(gelID) %>% unique() %>% pull()


l_pid_merged_main = unique(c(l_pid_ha_main, l_pid_eod_main, l_pid_hsp_main, l_pid_cmt_main, l_pid_eode_main, l_pid_als_main, l_pid_park_main, l_pid_cpark_main))
length(l_pid_merged_main)
# 3231
l_pid_merged_pilot = unique(c(l_pid_ha_pilot, l_pid_eod_pilot, l_pid_hsp_pilot, l_pid_cmt_pilot, l_pid_eode_pilot, l_pid_als_pilot, l_pid_park_pilot, l_pid_cpark_pilot))
length(l_pid_merged_pilot)
# 378

all_merged = unique(c(l_pid_merged_main,
                      l_pid_merged_pilot))
length(all_merged)
# 3609

# age
denak_main = table_a %>% filter(participant_id %in% all_merged) %>% select(age) %>% pull() 
denak_pilot = table_a_pilot %>% filter(gelID %in% all_merged) %>% select(age) %>% pull()
denak_age = c(denak_main, denak_pilot)

mean(denak_age)
# 53.5
summary(denak_age)
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#5.00   41.00   56.00   53.52   68.00  101.00 

# Gender
# Main
female_main = table_a %>% filter(participant_id %in% all_merged, participant_phenotypic_sex %in% "Female") %>% select(participant_id) %>% unique() %>% pull() %>% length()
male_main = table_a %>% filter(participant_id %in% all_merged, participant_phenotypic_sex %in% "Male") %>% select(participant_id) %>% unique() %>% pull() %>% length()

# Pilot
female_pilot = table_a_pilot %>% filter(gelID %in% all_merged, sex %in% "female") %>% select(gelID) %>% unique() %>% pull() %>% length()
male_pilot = table_a_pilot %>% filter(gelID %in% all_merged, sex %in% "male") %>% select(gelID) %>% unique() %>% pull() %>% length()

# Female (main + pilot) vs total
print((female_main + female_pilot) / length(all_merged))
# 0.468
# Male (main + pilot) vs total
print((male_main + male_pilot) / length(all_merged))
# 0.531

# Ethnicity
# MAIN
l_eth_main = table_a %>% filter(participant_id %in% all_merged) %>% select(participant_ethnic_category) %>% pull() 
# PILOT - consider all them `Not stated`
l_eth_pilot = rep("Not Stated", table_a_pilot %>% filter(gelID %in% all_merged) %>% select(gelID) %>% pull() %>% length())

l_eth_merged = c(l_eth_main, l_eth_pilot)
print(prop.table(table(l_eth_merged)))

####
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
# 112 25

# Simplify output TableA
table_a_expanded = table_a_expanded %>%
  select(list_samples, gene, allele, Repeat_Motif, participant_id, programme, genome_build, programme_consent_status, rare_diseases_family_id, biological_relationship_to_proband, 
         affection_status, participant_phenotypic_sex, year_of_birth, normalised_specific_disease, disease_sub_group, disease_group, family_group_type, family_medical_review_qc_state_code, 
         panel_list, best_guess_predicted_ancstry, self_reported, participant_ethnic_category, age, adult.paediatric)
colnames(table_a_expanded)[1] = "platekey" 
colnames(table_a_expanded)[3] = "repeat_size" 


write.table(table_a_expanded, "subtables/TableA_main.tsv", quote = F, row.names = F, col.names = T, sep = "\t")

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


write.table(table_a_pilot_expanded, "subtables/TableA_pilot.tsv", quote = F, row.names = F, col.names = T, sep = "\t")


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
# We need to take the list of PIDs from `list_2459_PIDs_ID_and_others_as_panels.txt`

# load list of 2449 PIDs 
l_complex_ID_group2 = read.table("list_2576_PIDs_ID_and_others_as_panels.txt", stringsAsFactors = F)
l_complex_ID_group2 = l_complex_ID_group2$V1
length(l_complex_ID_group2)
# 2576

table_b = table_diseases %>%
  filter(participant_id %in% l_complex_ID_group2)
dim(table_b)
# 2962  21

length(unique(table_b$plate_key.x))
# 2576


# Number of participants
l_pid_tableB = unique(table_b$participant_id)
length(l_pid_tableB)
# 2576

# Age
table_b %>% filter(participant_id %in% l_pid_tableB) %>% select(age) %>% pull %>% mean()
# 15.03
table_b %>% filter(participant_id %in% l_pid_tableB) %>% select(age) %>% pull %>% summary()
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 3.00    8.00   12.00   15.03   19.00   90.00 

# Females
table_b %>% filter(participant_id %in% l_pid_tableB, participant_phenotypic_sex %in% "Female") %>% select(participant_id) %>% unique() %>% pull() %>% length()
# 1141
table_b %>% filter(participant_id %in% l_pid_tableB, participant_phenotypic_sex %in% "Male") %>% select(participant_id) %>% unique() %>% pull() %>% length()
# 1435

# Ethnicity
table_b$participant_ethnic_category = recode(table_b$participant_ethnic_category,
                                             "White: British"= "White",
                                             "White White: British"= "White",
                                             "White: Any other White background"="White",
                                             "Asian or Asian British: Pakistani"="Asian",
                                             "Asian or Asian British: Indian"="Asian",
                                             "Asian or Asian British: Any other Asian background"="Asian",
                                             "Black or Black British: African"="Black",
                                             "Other Ethnic Groups: Any other ethnic group"="Other", 
                                             "Mixed: Any other mixed background"="Mixed",
                                             "White: Irish"="White",
                                             "Asian or Asian British: Bangladeshi"="Asian",
                                             "Mixed: White and Asian"="Mixed",
                                             "Mixed: White and Black Caribbean"="Mixed",
                                             "Black or Black British: Caribbean"="Black",
                                             "Mixed: White and Black African"="Mixed",
                                             "Black or Black British: Any other Black background"="Black",
                                             "Other Ethnic Groups: Chinese"="Other")
# Defining NA's as `Not stated`
which_na = which(is.na(table_b$participant_ethnic_category))
table_b$participant_ethnic_category[which_na] = "Not Stated"

table_b = unique(table_b)

table_b %>% filter(participant_id %in% l_pid_tableB) %>% select(participant_ethnic_category) %>% table() %>% prop.table()
#     Asian       Black       Mixed  Not Stated       Other       White 
# 0.131667792 0.022282242 0.034098582 0.169480081 0.009453072 0.633018231 

# Let's define the list of genes for Table B
l_genes_tableB = c("ATN1_CAG","ATXN1_CAG", "ATXN2_CAG", "ATXN3_CAG", "ATXN7_CAG", "HTT_CAG", "TBP_CAG")

l_platekeys_tableB = unique(table_b$plate_key.x)

# Define specific thresholds for TABLE B
gene_pathogenic_threshold_tableB = data.frame(locus = l_genes_tableB,
                                              threshold = c(63,44,60,75,90,60,60))

# Now, we want to see how many of them have an expansion on any of the genes in B
expanded_table_main = data.frame()
for (i in 1:length(l_genes_tableB)){
  locus_name = l_genes_tableB[i]
  patho_cutoff = gene_pathogenic_threshold_tableB %>% 
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
# 28  5

# separate platekeys
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
# 187  5

# From the expanded table, let's see how many are in l_platekeys_tableB
expanded_table_main_in_tableB = expanded_table_main_per_locus %>%
  filter(list_samples %in% l_platekeys_tableB)
dim(expanded_table_main_in_tableB)
# 11  5

# Let' enrich expanded TABLE B repeats with clinical data from `table_b`
table_b_expanded = left_join(expanded_table_main_in_tableB,
                             table_b,
                             by = c("list_samples" = "plate_key.x"))
dim(table_b_expanded)
# 12  25

# Simplify output TableB
table_b_expanded = table_b_expanded %>%
  select(list_samples, gene, allele, Repeat_Motif, participant_id, programme, genome_build, programme_consent_status, rare_diseases_family_id, biological_relationship_to_proband, 
         affection_status, participant_phenotypic_sex, year_of_birth, normalised_specific_disease, disease_sub_group, disease_group, family_group_type, family_medical_review_qc_state_code, 
         panel_list, best_guess_predicted_ancstry, self_reported, participant_ethnic_category, age, adult.paediatric)
colnames(table_b_expanded)[1] = "platekey" 
colnames(table_b_expanded)[3] = "repeat_size" 
write.table(table_b_expanded, "subtables/TableB_main.tsv", quote = F, row.names = F, col.names = T, sep = "\t")

# This is the raw data for Table B - Main
# We want to compute all diseases in the coctail of `l_diseases_tableB` as 1
l_diseases_tableB = unique(table_b$normalised_specific_disease)
matrix_to_print = matrix(nrow  = length(l_diseases_tableB), ncol = length(l_genes_tableB))
for(i in 1:length(l_diseases_tableB)){
  for (j in 1:length(l_genes_tableB)){
    number_to_print = table_b_expanded %>% 
      filter(normalised_specific_disease %in% l_diseases_tableB[i], gene %in% l_genes_tableB[j]) %>% 
      select(participant_id) %>% unique() %>% pull() %>% length()
    
    print(l_diseases_tableB[i])
    print(l_genes_tableB[j])
    print(number_to_print)
    matrix_to_print[i,j] = number_to_print
  }
}

rownames(matrix_to_print) = l_diseases_tableB
colnames(matrix_to_print) = l_genes_tableB

write.table(matrix_to_print, "./subtables/tableB_main_for_excel.tsv", sep = "\t", row.names = T, col.names = T, quote = F)


################################################################################################################################################################
# TABLE C
# patients presenting with intellectual disability and or a neuromuscular phenotype were analysed for DMPK

# As for DMPK: select patients in this way: 
# 1) recruited under intellectual disability (normalised spec disease) AND who have been applied EITHER one of the following panel: 
# "Congenital muscular dystrophy", " Congenital myopathy", "Skeletal Muscle Channelopathies"; 

# 2) adult and children only that are recruited under specific disease "Congenital muscular dystrophy" OR 
# " Congenital myopathy" OR  
# "Skeletal Muscle Channelopathies" OR 
# "Distal myopathies": 

# Function that checks if any of the items in list of characters 1 does exist in list of characters 2
any_exist <- function(list1, list2) {
  for (i in list1){
    if (i %in% list2){
      return(TRUE)
    }
  }
  return(FALSE)
}

# Since we are going to work with panels, let's make life simple
table_panels_row = table_diseases %>% 
  select(plate_key.x, participant_id, normalised_specific_disease, disease_sub_group, disease_group, panel_list) %>%
  mutate(panels = strsplit(as.character(panel_list), ",")) %>%
  unnest(panels) %>%
  as.data.frame()
table_panels_row = unique(table_panels_row)
dim(table_panels_row)
# 52302  7

table_panels_row_pilot = table_diseases_pilot %>%
  select(plateKey, gelID, specificDisease, panel_list) %>%
  mutate(panels = strsplit(as.character(panel_list), ",")) %>%
  unnest(panels) %>%
  as.data.frame()
table_panels_row_pilot = unique(table_panels_row_pilot)
dim(table_panels_row_pilot)
# 1129  5


# We will create 2 groups as above
l_group1 = c("Congenital muscular dystrophy", "")
# MAIN
table_c = table_diseases %>%
  filter(normalised_specific_disease %in% c("Intellectual disability",
                                            "Kabuki syndrome",
                                            "Congenital muscular dystrophy",
                                            "Congenital myopathy",
                                            "Skeletal Muscle Channelopathies",
                                            "Distal myopathies"))
dim(table_c)
# 7695  21

# PILOT
table_c_pilot = table_diseases_pilot %>%
  filter(specificDisease %in% c("Intellectual disability",
                                "Kabuki syndrome",
                                "Congenital muscular dystrophy",
                                "Congenital myopathy",
                                "Skeletal Muscle Channelopathies",
                                "Distal myopathies"))
dim(table_c_pilot)
# 242  15

# Let's define the list of genes for Table C
l_genes_tableC = c("DMPK_CTG")

# How many PIDs in the Main?
length(unique(table_c$participant_id))
# 7345

# How many PIDs are in the Pilot?
length(unique(table_c_pilot$plateKey))
# 241

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

# separate platekeys
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
write.table(table_c_expanded, "subtables/TableC_main.tsv", quote = F, row.names = F, col.names = T, sep = "\t")

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
# TABLE D. only including children recruited under ID (using >55. as cutoff)	
# FMR1
# Intellectual disability	

# MAIN
table_d = table_diseases %>%
  filter(normalised_specific_disease %in% c("Intellectual disability",
                                            "Kabuki syndrome"))
dim(table_d)
# 6890  21

# PILOT
table_d_pilot = table_diseases_pilot %>%
  filter(specificDisease %in% c("Intellectual disability",
                                "Kabuki syndrome"))
dim(table_d_pilot)
# 161  15

# Recode ethnicity
table_d$participant_ethnic_category = recode(table_d$participant_ethnic_category,
                                             "White: British"= "White",
                                             "White White: British"= "White",
                                             "White: Any other White background"="White",
                                             "Asian or Asian British: Pakistani"="Asian",
                                             "Asian or Asian British: Indian"="Asian",
                                             "Asian or Asian British: Any other Asian background"="Asian",
                                             "Black or Black British: African"="Black",
                                             "Other Ethnic Groups: Any other ethnic group"="Other", 
                                             "Mixed: Any other mixed background"="Mixed",
                                             "White: Irish"="White",
                                             "Asian or Asian British: Bangladeshi"="Asian",
                                             "Mixed: White and Asian"="Mixed",
                                             "Mixed: White and Black Caribbean"="Mixed",
                                             "Black or Black British: Caribbean"="Black",
                                             "Mixed: White and Black African"="Mixed",
                                             "Black or Black British: Any other Black background"="Black",
                                             "Other Ethnic Groups: Chinese"="Other")
# Defining NA's as `Not stated`
which_na = which(is.na(table_d$participant_ethnic_category))
table_d$participant_ethnic_category[which_na] = "Not Stated"

table_d = unique(table_d)
table_d_pilot = unique(table_d_pilot)

# Num participants
length(unique(table_d$participant_id))
# 6570
length(unique(table_d_pilot$gelID))
# 161

# Age
l_age_main = table_d$age
l_age_pilot = table_d_pilot$age
l_age_merged = c(l_age_main, l_age_pilot)

mean(l_age_merged)
# 13.5
summary(l_age_merged)
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#2.00    8.00   11.00   13.49   16.00   72.00 


# GEnder
table_d %>% filter(participant_phenotypic_sex %in% "Female") %>% select(participant_id) %>% unique() %>% pull() %>% length()
# 2615
table_d %>% filter(participant_phenotypic_sex %in% "Male") %>% select(participant_id) %>% unique() %>% pull() %>% length()
# 3955

table_d_pilot %>% filter(sex %in% "female") %>% select(gelID) %>% unique() %>% pull() %>% length()
# 65
table_d_pilot %>% filter(sex %in% "male") %>% select(gelID) %>% unique() %>% pull() %>% length()
# 96

(2615 + 65) / 6731
# 0.395
(3955 + 96) / 6731
# 0.60

# Ethnicity
l_eth_main = table_d %>% select(participant_ethnic_category) %>% pull() 
# PILOT - consider all them `Not stated`
l_eth_pilot = rep("Not Stated", length(unique(table_d_pilot$gelID)))

l_eth_merged = c(l_eth_main, l_eth_pilot)
print(prop.table(table(l_eth_merged)))
# Asian      Black      Mixed Not Stated      Other      White
#0.08920720 0.01942987 0.03616508 0.19472415 0.01049497 0.64997873 

###
# Let's define the list of genes for Table C
l_genes_tableD = c("FMR1_CGG")

# How many PIDs in the Main?
length(unique(table_d$participant_id))
# 6570

# How many PIDs are in the Pilot?
length(unique(table_d_pilot$plateKey))
# 161

l_platekeys_tableD = unique(table_d$plate_key.x)
l_platekeys_tableD_pilot = unique(table_d_pilot$plateKey)

# Now, we want to see how many of them have an expansion on any of the genes in `FMR1` (cutoff >= 55)
expanded_table_main = data.frame()
for (i in 1:length(l_genes_tableD)){
  locus_name = l_genes_tableD[i]
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
# 113  5

# Now, we want to see how many of them have an expansion on any of the genes in `l_genes_tableA` - but for Pilot data
expanded_table_pilot = data.frame()
for (i in 1:length(l_genes_tableD)){
  locus_name = l_genes_tableD[i]
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
# 21  5

# separate platekeys
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
# 2147  5

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

# From the expanded table, let's see how many are in l_platekeys_tableC
expanded_table_main_in_tableD = expanded_table_main_per_locus %>%
  filter(list_samples %in% l_platekeys_tableD)
dim(expanded_table_main_in_tableD)
# 159  5

# The same por PILOT
expanded_table_pilot_in_tableD = expanded_table_pilot_per_locus %>%
  filter(list_samples %in% l_platekeys_tableD_pilot)
dim(expanded_table_pilot_in_tableD)
# 0  5

# Let' enrich expanded TABLE C repeats with clinical data from `table_a`
table_d_expanded = left_join(expanded_table_main_in_tableD,
                             table_d,
                             by = c("list_samples" = "plate_key.x"))
dim(table_d_expanded)
# 170  25

# PILOT - nothing to merge

# Simplify output TableD
table_d_expanded = table_d_expanded %>%
  select(list_samples, gene, allele, Repeat_Motif, participant_id, programme, genome_build, programme_consent_status, rare_diseases_family_id, biological_relationship_to_proband, 
         affection_status, participant_phenotypic_sex, year_of_birth, normalised_specific_disease, disease_sub_group, disease_group, family_group_type, family_medical_review_qc_state_code, 
         panel_list, best_guess_predicted_ancstry, self_reported, participant_ethnic_category, age, adult.paediatric)
colnames(table_d_expanded)[1] = "platekey" 
colnames(table_d_expanded)[3] = "repeat_size" 

dim(table_d_expanded)
# 170  24

# we need to only include children recruited under ID
table_d_expanded = table_d_expanded %>%
  filter(adult.paediatric %in% "Paediatric")
dim(table_d_expanded)
# 132 24

write.table(table_d_expanded, "subtables/TableD_main.tsv", quote = F, row.names = F, col.names = T, sep = "\t")

# This is the raw data for Table D- Main
# Let's do numbers for FMR1 and Intellectual disability
l_diseases_tableD = unique(table_d$normalised_specific_disease)
matrix_to_print = matrix(nrow  = length(l_diseases_tableD), ncol = 1)
for(i in 1:length(l_diseases_tableD)){
  for (j in 1:length(l_genes_tableD)){
    number_to_print = table_d_expanded %>% 
      filter(normalised_specific_disease %in% l_diseases_tableD[i], gene %in% l_genes_tableD[j]) %>% 
      select(participant_id) %>% unique() %>% pull() %>% length()
    
    print(l_diseases_tableD[i])
    print(l_genes_tableD[j])
    print(number_to_print)
    matrix_to_print[i,j] = number_to_print
  }
}

rownames(matrix_to_print) = l_diseases_tableD
colnames(matrix_to_print) = l_genes_tableD
write.table(matrix_to_print, "./subtables/tableD_main_for_excel.tsv", sep = "\t", row.names = T, col.names = T, quote = F)







