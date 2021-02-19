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
setwd("~/Documents/STRs/ANALYSIS/population_research/PAPER/expanded_genomes_main_pilot/feb2021/")

# load data
merged_table = read.csv("~/Documents/STRs/data/research/batch_august2020/output_EHv3.2.2_vcfs/merged/merged_93446_genomes_EHv322_batch_august2020.tsv",
                        sep = "\t",
                        stringsAsFactors = F, 
                        header = T)
dim(merged_table)
# 27238  12

# Focus on our 13 loci
merged_table = merged_table %>%
  filter(gene %in% c("AR", "ATN1","ATXN1", "ATXN2", "ATXN3", "ATXN7", "C9ORF72",
                     "CACNA1A", "DMPK", "FMR1", "FXN", "HTT", "TBP"))

# Load clinical data (Main and Pilot) - even though we will keep with unrel genomes
clin_data = read.csv("~/Documents/STRs/clinical_data/clinical_data/Main_RE_V10_and_Pilot_programmes.tsv",
                      stringsAsFactors = F, 
                      sep = "\t",
                      header = T)
dim(clin_data)
# 2101385  24

# Load unrelated list of genomes 
l_unrelated = read.table("~/Documents/STRs/ANALYSIS/population_research/MAIN_ANCESTRY/batch2/l_unrelated_55603_genomes_batch2.txt",
                         stringsAsFactors = F)
l_unrelated = l_unrelated$V1
length(l_unrelated)
# 55603

# Disease_group info (and all other clinical characteristics) we've got for probands
# Let's take the familyIDs that have been recruited as Neuro in `disease_group`
l_fam_neuro = clin_data %>%
  filter(grepl("[Nn][Ee][Uu][Rr][Oo]", diseasegroup_list)) %>%
  select(rare_diseases_family_id) %>%
  unique() %>%
  pull()
length(l_fam_neuro)
# 14402

clin_data = clin_data %>% select(platekey, rare_diseases_family_id, diseasegroup_list, superpopu)
clin_data = unique(clin_data)
dim(clin_data)
# 109411  4

# Load platekey-pid-famID table we created to fish platekeys not included in further RE releases
clin_metadata = read.csv("~/Documents/STRs/clinical_data/clinical_data/merged_RE_releases_and_Pilot_PID_FID_platekey.tsv",
                         stringsAsFactors = F,
                         sep = "\t",
                         header = T)
dim(clin_metadata)
# 621704  4

# Include or enrich `clin_data` with extra platekeys, to associate platekey <-> famID
clin_data = full_join(clin_data,
                      clin_metadata %>% select(platekey, participant_id, rare_diseases_family_id),
                      by = "platekey")
clin_data = unique(clin_data)
dim(clin_data)
# 149776  6

# First let's unite `rare_diseases_family_id` columns into 1
clin_data = clin_data %>%
  group_by(rare_diseases_family_id.x) %>%
  mutate(famID = ifelse(is.na(rare_diseases_family_id.x), rare_diseases_family_id.y, rare_diseases_family_id.x)) %>%
  ungroup() %>%
  as.data.frame()

# Now we've got complete famID, let's define whether each platkey is neuro or not
clin_data = clin_data %>% 
  group_by(famID) %>% 
  mutate(is_neuro = ifelse(famID %in% l_fam_neuro, "Neuro", "NotNeuro")) %>% 
  ungroup() %>% 
  as.data.frame() %>%
  select(platekey, participant_id, famID, diseasegroup_list, is_neuro)

# Let's include a column which says whether a platekey is unrel or not
clin_data = clin_data %>%
  group_by(platekey) %>%
  mutate(is_unrel = ifelse(platekey %in% l_unrelated, TRUE, FALSE)) %>%
  ungroup() %>%
  as.data.frame()

# Check if we have 55,603 unrel genomes
clin_data %>% filter(is_unrel) %>% select(platekey) %>% unique() %>% pull() %>% length()
# 55603

# Check whether all 55,603 unrel genomes we want to analyse and will consider as TOTAL DATASET, are also included in the EHv322 table
l_all_samples_merged = c()
for (i in 1:length(merged_table$chr)){
  aux_vcf = strsplit(merged_table$list_samples[i], ';')[[1]]
  aux_vcf = gsub("EH_", "", aux_vcf)
  aux_vcf = gsub(".vcf", "", aux_vcf)
  aux_vcf = gsub("_x2", "", aux_vcf)
  l_all_samples_merged = unique(c(l_all_samples_merged,
                           aux_vcf))
}
length(l_all_samples_merged)
# 93446

# How many intersect, difference?
length(intersect(l_all_samples_merged, l_unrelated))
# 55,599
length(setdiff(l_unrelated, l_all_samples_merged))
# 4
# "LP3000170-DNA_D10" "LP3000170-DNA_D11" "LP3000170-DNA_D12" "LP3000115-DNA_G11"

#write down  these extra 644 in order to include in the next run
write.table(setdiff(l_unrelated_merged, l_all_samples_merged),
            "~/Documents/STRs/data/research/input/list_644_genomes_not_included_in_batch_march2020_but_popu.tsv",
            quote = F,
            row.names = F,
            col.names = F)

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


# Take the ones that intersect and see how many genomes from each superpopu we have
l_genomes_batch_and_popu = intersect(l_all_samples_merged, l_unrelated_merged)
popu_table %>% filter(ID %in% l_genomes_batch_and_popu) %>% select(merged_superpopu) %>% pull() %>% table()
#  AFR   AMR   EAS   EUR   SAS 
# 1599  1074   177 32207  3361 

pilot_popu_table %>% filter(ID %in% l_genomes_batch_and_popu) %>% select(merged_superpopu_pilot) %>% pull() %>% table()
# AFR  AMR  EAS  EUR  SAS 
# 43   25    4 1779  161 

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
    filter(ID %in% l_unrelated_merged)

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





# Beyond full-mutation or pathogenic cutoff
l_genes = c("AR", "ATN1", "ATXN1", "ATXN2", "ATXN3", "ATXN7", "C9ORF72", "CACNA1A", "DMPK", "FMR1", "FXN", "HTT", "TBP")
l_patho_cutoff = c(38,48,44,33,60,36,60,20,50,200,66,40,49)
l_list_unrelated = list.files(pattern = "_unrelated.tsv$")
for (i in 1:length(l_genes)){
  print(l_genes[i])
  locus_table = read.csv(l_list_unrelated[i],
                         stringsAsFactors = F,
                         header = T,
                         sep = "\t")
  locus_table %>% filter(repeat_size >= l_patho_cutoff[i]) %>% select(ID) %>% unique() %>% pull() %>% length() %>% print()
}

