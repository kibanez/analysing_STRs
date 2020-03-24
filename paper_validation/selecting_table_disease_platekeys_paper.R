# Objective: select the platekeys/participant ids corresponding to the list of diseases we consider for the paper
# libraries
library(dplyr)
library(tidyverse); packageDescription ("tidyverse", fields = "Version") # "1.2.1

setwd("~/Documents/STRs/PAPER/VALIDATION/")

# `clinical_data_research_cohort_86457_genomes_withPanels_250919` is a table I generated

table_diseases <- clinical_data_research_cohort_86457_genomes_withPanels_250919 %>%
clinical_data_research_cohort_86457_genomes_withPanels_250919 %>% count(normalised_specific_disease)

count_of_diseases <- clinical_data_research_cohort_86457_genomes_withPanels_250919 %>% count(normalised_specific_disease)

table_diseases <- clinical_data_research_cohort_86457_genomes_withPanels_250919 %>% 
  filter(normalised_specific_disease %in% c("Intellectual disability",
                                            "Amyotrophic lateral sclerosis or motor neuron disease", 
                                            "Charcot-Marie-Tooth disease", "Congenital muscular dystrophy",
                                            "Congenital myopathy", "Early onset dementia", "Early onset dystonia", 
                                            "Distal myopathies", "Complex Parkinsonism", "Hereditary ataxia", 
                                            "Hereditary spastic paraplegia", 
                                            "'Early onset and familial Parkinson''s Disease'"))
dim(table_diseases)
table(table_diseases$normalised_specific_disease)

parkinson_to_enrich = clinical_data_research_cohort_86457_genomes_withPanels_250919 %>% 
  filter(grepl("Complex Parkin", normalised_specific_disease))
dim(parkinson_to_enrich)
table(parkinson_to_enrich$normalised_specific_disease)

table_diseases = rbind(table_diseases, parkinson_to_enrich)
dim(table_diseases)

write.table(table_diseases[ ,2], row.names = FALSE, col.names = FALSE, quote = FALSE, "~/Documents/STR identif/clinical data lp numbers/list_PIDs_table_diseases.txt")


pid_platekey_genomeBuild_table_diseases_dedup <- read.csv("~/Documents/STR identif/clinical data lp numbers/pid_platekey_genomeBuild_table_diseases_dedup.csv", stringsAsFactors=FALSE)

list_197_PIDs_from_table_diseases <- read.table("~/Documents/STR identif/clinical data lp numbers/list_197_PIDs_from_table_diseases.txt", quote="\"", comment.char="", stringsAsFactors=FALSE)


list_197_PIDs_from_table_diseases = list_197_PIDs_from_table_diseases$V1

list_platekeys = clinical_data_research_cohort_86457_genomes_withPanels_250919 %>% 
  filter(participant_id %in% list_197_PIDs_from_table_diseases) %>% select(participant_id, plate_key.x, genome_build)
dim(list_platekeys)

list_platekeys = list_platekeys %>% group_by(participant_id) %>% mutate(latest = max(plate_key.x)) %>% ungroup() %>% as.data.frame()

list_platekeys = list_platekeys %>% select(participant_id, latest, genome_build)
list_platekeys = unique(list_platekeys)
dim(list_platekeys)

l_platekey_dedup = list_platekeys$participant_id[which(duplicated(list_platekeys$participant_id))]
length(l_platekey_dedup)

list_platekeys = list_platekeys %>% filter(!(participant_id %in% l_platekey_dedup & genome_build %in% "GRCh37"))
dim(list_platekeys)

colnames(list_platekeys) = colnames(pid_platekey_genomeBuild_table_diseases_dedup)

pid_platekey_genomeBuild_table_diseases_dedup = rbind(pid_platekey_genomeBuild_table_diseases_dedup, list_platekeys)
dim(pid_platekey_genomeBuild_table_diseases_dedup)

table_diseases_enriched = clinical_data_research_cohort_86457_genomes_withPanels_250919 %>% 
  filter(plate_key.x %in% pid_platekey_genomeBuild_table_diseases_dedup$platekey)
dim(table_diseases_enriched)

table_diseases_enriched = unique(table_diseases_enriched)
dim(table_diseases_enriched)

write.table(table_diseases_enriched, file = "~/Documents/STR identif/clinical data lp numbers/table_diseases_enriched.tsv", sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)


