# Check which participants we are taking into account in the population research
# are they all unrelated, probands, affected one?
date ()
Sys.info ()[c("nodename", "user")]
commandArgs ()
rm (list = ls ())
R.version.string ## "R version 3.3.2 (2016-10-31)"

# libraries
library(dplyr); packageDescription ("dplyr", fields = "Version") #"0.8.3"
library(ggplot2); packageDescription ("ggplot2", fields = "Version") #"3.2.1"

# set the working directory
setwd("/Users/kibanez/Documents/STRs/ANALYSIS/population_research")

l_participants_pure_ancestry=read.table("list_genomes_56176_pure_ancestry_info.txt",
                                      header = F,
                                      stringsAsFactors = F)
l_participants_pure_ancestry = l_participants_pure_ancestry$V1
length(l_participants_pure_ancestry)
# 56176


re_data = read.csv("~/Documents/STRs/clinical_data/clinical_data/rd_genomes_all_data_250919.tsv",
                   header = T,
                   stringsAsFactors = F,
                   sep = "\t")
dim(re_data)
# 1056568  26

re_data_subset = re_data %>%
  filter(plate_key.x %in% l_participants_pure_ancestry) %>%
  select(participant_id, plate_key.x, rare_diseases_family_id, participant_type, normalised_specific_disease, genome_build, disease_group, disease_sub_group, specific_disease, programme, family_group_type, affection_status)
dim(re_data_subset)
# 634692  12

length(unique(re_data_subset$plate_key.x))
# 56176

dim(unique(re_data_subset))
# 57704

re_data_subset = unique(re_data_subset)

table(re_data_subset$programme)
# Cancer Rare Diseases 
# 8655         49049 

# We want unrelated individuals - probands
re_data_subset_unrelatd_probands = re_data_subset %>%
  filter((programme %in% "Cancer") | (programme %in% "Rare Diseases" & participant_type %in% "Proband")) 
dim(re_data_subset_unrelatd_probands)  
# 32440  12

length(unique(re_data_subset_unrelatd_probands$plate_key.x))
# 31336

# Curiosity -  how many have neuro as `disease_group`
table(re_data_subset_unrelatd_probands$disease_group)

# Let's enrich them with population ancestry information
l_AFR = read.table("./list_genomes_1794_pure_AFR.txt", stringsAsFactors = F)
l_AMR = read.table("./list_genomes_803_pure_AMR.txt", stringsAsFactors = F)
l_ASI = read.table("./list_genomes_5971_pure_ASI.txt", stringsAsFactors = F)
l_EAS = read.table("./list_genomes_400_pure_EAS.txt", stringsAsFactors = F)
l_EUR = read.table("./list_genomes_47208_pure_EUR.txt", stringsAsFactors = F)

l_AFR = l_AFR$V1
length(l_AFR)
# 1794

l_AMR = l_AMR$V1
length(l_AMR)
# 803

l_ASI = l_ASI$V1
length(l_ASI)
# 5971

l_EAS = l_EAS$V1
length(l_EAS)
# 400

l_EUR = l_EUR$V1
length(l_EUR)
# 47208

l_AFR_unrelated_probands = re_data_subset_unrelatd_probands %>%
  filter(plate_key.x %in% l_AFR) %>%
  select(plate_key.x) %>% unique() %>% pull() %>% as.character()
length(l_AFR_unrelated_probands)
# 1164

l_AMR_unrelated_probands = re_data_subset_unrelatd_probands %>%
  filter(plate_key.x %in% l_AMR) %>%
  select(plate_key.x) %>% unique() %>% pull() %>% as.character()
length(l_AMR_unrelated_probands)
# 396

l_ASI_unrelated_probands = re_data_subset_unrelatd_probands %>%
  filter(plate_key.x %in% l_ASI) %>%
  select(plate_key.x) %>% unique() %>% pull() %>% as.character()
length(l_ASI_unrelated_probands)
# 2761

l_EAS_unrelated_probands = re_data_subset_unrelatd_probands %>%
  filter(plate_key.x %in% l_EAS) %>%
  select(plate_key.x) %>% unique() %>% pull() %>% as.character()
length(l_EAS_unrelated_probands)
# 231

l_EUR_unrelated_probands = re_data_subset_unrelatd_probands %>%
  filter(plate_key.x %in% l_EUR) %>%
  select(plate_key.x) %>% unique() %>% pull() %>% as.character()
length(l_EUR_unrelated_probands)
# 26784

# Let's write down the list of platekeys that correspond to unrelated, proband OR cancer genomes, for each ancestry
write.table(l_AFR_unrelated_probands, "./list_genomes_unrelated_proband_or_cancer_1164_pure_AFR.txt", quote = F, row.names = F, col.names = F)
write.table(l_AMR_unrelated_probands, "./list_genomes_unrelated_proband_or_cancer_396_pure_AMR.txt", quote = F, row.names = F, col.names = F)
write.table(l_ASI_unrelated_probands, "./list_genomes_unrelated_proband_or_cancer_2761_pure_ASI.txt", quote = F, row.names = F, col.names = F)
write.table(l_EAS_unrelated_probands, "./list_genomes_unrelated_proband_or_cancer_231_pure_EAS.txt", quote = F, row.names = F, col.names = F)
write.table(l_EUR_unrelated_probands, "./list_genomes_unrelated_proband_or_cancer_26784_pure_EUR.txt", quote = F, row.names = F, col.names = F)


