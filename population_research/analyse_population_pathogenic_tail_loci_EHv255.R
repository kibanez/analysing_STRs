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
setwd("~/Documents/STRs/ANALYSIS/population_research/PAPER/expanded_genomes_main_pilot/feb2021/EHv255/")

# load data
merged_table = read.csv("~/Documents/STRs/data/research/batch_august2020/output_EHv2.5.5_vcfs/merged/merged_93446_genomes_EHv255_batch_august2020.tsv",
                        sep = "\t",
                        stringsAsFactors = F, 
                        header = T)
dim(merged_table)
# 32703  11

# Focus on our 13 loci
merged_table = merged_table %>%
  filter(gene %in% c("AR_CAG", "ATN1_CAG","ATXN1_CAG", "ATXN2_CAG", "ATXN3_CAG", "ATXN7_CAG", "C9orf72_GGGGCC",
                     "CACNA1A_CAG", "DMPK_CAG", "FMR1_CGG", "FXN_GAA", "HTT_CAG", "TBP_CAG"))

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
# 122024  4

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
  select(platekey, participant_id, famID, diseasegroup_list, is_neuro, superpopu)

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
  aux_vcf = gsub("_x2", "", aux_vcf)
  aux_vcf = gsub(".vcf", "", aux_vcf)
  
  l_all_samples_merged = unique(c(l_all_samples_merged,
                           aux_vcf))
}
length(l_all_samples_merged)
# 93535

# How many intersect, difference?
length(intersect(l_all_samples_merged, l_unrelated))
# 55,171
length(setdiff(l_unrelated, l_all_samples_merged))
# 432 - I've realised all these are empty!!! I need to re-run EHv255 on them.

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

# For each locus
l_genes = c("AR", "ATN1", "ATXN1", "ATXN2", "ATXN3", "ATXN7","C9ORF72", "CACNA1A", "DMPK", "FMR1", "FXN", "HTT", "TBP")
l_premut_cutoff = c(34,34,39,31,43,34,30,17,50,55,44,35,41)
l_patho_cutoff = c(38,48,44,33,60,36,60,20,50,200,66,40,49)

df_cutoff = data.frame(locus = l_genes,
                       premut_cutoff = l_premut_cutoff,
                       patho_cutoff = l_patho_cutoff)

for (i in 1:length(l_genes)){
  print(l_genes[i])
  
  # Beyond premut
  merged_table_locus = merged_table %>%
    filter(gene %in% l_genes[i], allele > l_premut_cutoff[i])

  list_allele_size = c()
  list_vcf_premut_locus = c()
  df_platekey_size = data.frame()
  for (j in 1:length(merged_table_locus$list_samples)){
    list_vcf_allele = strsplit(merged_table_locus$list_samples[j], ';')[[1]]
    number_vcf = length(list_vcf_allele)
    list_vcf_premut_locus = c(list_vcf_premut_locus,
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
  
  list_vcf_premut_locus = gsub('.vcf', '', list_vcf_premut_locus)
  list_vcf_premut_locus = gsub('^EH_', '', list_vcf_premut_locus)
  list_vcf_premut_locus = gsub('_x2', '', list_vcf_premut_locus)
  
  # Beyond patho threshold
  merged_table_locus = merged_table %>%
    filter(gene %in% l_genes[i], allele >= l_patho_cutoff[i])
  
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
  
  # Enrich platekeys with simplified clinical data: popu, is_unrel, neuro_notNeuro
  premut_popu = clin_data %>%
    filter(platekey %in% list_vcf_premut_locus)
  premut_popu = unique(premut_popu)
  print(dim(premut_popu))

  patho_popu = clin_data %>%
    filter(platekey %in% list_vcf_patho_locus)
  patho_popu = unique(patho_popu)
  print(dim(patho_popu))
  
  # Add locus name as column
  premut_popu$locus = rep(l_genes[i], length(premut_popu$platekey))
  patho_popu$locus = rep(l_genes[i], length(patho_popu$platekey))
  
  output_file_name = paste(l_genes[i], "beyond_", sep = "_")
  
  output_file_name1 = paste(output_file_name, "premutation_cutoff_", sep = "_")
  output_file_name1 = paste(output_file_name1, as.character(l_premut_cutoff[i]), sep = "")
  output_file_name1 = paste(output_file_name1, "EHv322_92K_population.tsv", sep = "_")
  output_file_name1 = paste("./beyond_premut/", output_file_name1, sep = "")
  
  output_file_name2 = paste(output_file_name, "pathogenic_cutoff_", sep = "_")
  output_file_name2 = paste(output_file_name2, as.character(l_patho_cutoff[i]), sep = "")
  output_file_name2 = paste(output_file_name2, "EHv322_92K_population.tsv", sep = "_")
  output_file_name2 = paste("./beyond_full-mutation/", output_file_name2, sep = "")
  
  write.table(premut_popu, 
              output_file_name1, 
              sep = "\t",
              quote = F,
              row.names = F,
              col.names = T)
  
  write.table(patho_popu, 
              output_file_name2, 
              sep = "\t",
              quote = F,
              row.names = F,
              col.names = T)
}
  
