# Objective: analyse haplotyping blocks within HTT
date ()
Sys.info ()[c("nodename", "user")]
commandArgs ()
rm (list = ls ())
R.version.string ## "R version 3.6.1 (2019-07-05)"

# libraries
library(dplyr)

# Set working directory
setwd("~/Documents/STRs/ANALYSIS/haplotyping/HTT/feb2021/gvcfgenotyper/")

# Load clinical data
clin_data = read.csv("~/Documents/STRs/clinical_data/clinical_data/rd_genomes_all_data_251120_V10.tsv",
                     sep = "\t",
                     stringsAsFactors = F,
                     header = T)
dim(clin_data)
# 2096500  36

# Load 29 expanded genomes after visual inspection and unrelated
l_exp_genomes = read.table("~/Documents/STRs/ANALYSIS/population_research/100K/carrier_freq/list_29_expanded_after_QC_unrelated.tsv", stringsAsFactors = F)
l_exp_genomes = l_exp_genomes$V1
length(l_exp_genomes)
# 29

# Retrieve genome build, sex, popu for these genomes
haplo_genomes = clin_data %>%
  filter(platekey %in% l_exp_genomes) %>%
  select(platekey, genome_build, participant_ethnic_category, participant_phenotypic_sex, superpopu, programme, file_path)
haplo_genomes = unique(haplo_genomes)
dim(haplo_genomes)
# 29  7

length(unique(haplo_genomes$platekey))
# 29

# Adapt or change path to the BAM file, to path to the gVCF file
haplo_genomes$file_path = gsub("Assembly","Variations", haplo_genomes$file_path)
haplo_genomes$file_path = gsub("bam","genome.vcf.gz", haplo_genomes$file_path)

# write them into a file
write.table(haplo_genomes,
            "./table_29genomes_expanded_unrelated_HTT.tsv",
            sep = "\t",
            quote = F,
            row.names = F,
            col.names = T)

# Also write here, the phenotype-part for CASES for later
# ID, IID, sex, CaseControl
pheno_cases = haplo_genomes %>% 
  select(platekey, participant_phenotypic_sex)
pheno_cases$ID = seq(1,length(pheno_cases$platekey), 1)
pheno_cases$ID = paste("FAM", pheno_cases$ID, sep = '_')
pheno_cases$CaseControl = rep("2", length(pheno_cases$platekey))

# format females as `2` and males as `1`
pheno_cases$participant_phenotypic_sex = gsub("Female", "2", pheno_cases$participant_phenotypic_sex)
pheno_cases$participant_phenotypic_sex = gsub("Male", "1", pheno_cases$participant_phenotypic_sex)

# Reorder columns and rename them
pheno_cases = pheno_cases %>%
  select(ID, platekey, participant_phenotypic_sex, CaseControl)
colnames(pheno_cases) = c("ID", "IID", "sex", "CaseControl")
write.table(pheno_cases, "list_phenotypes_cases.tsv", quote = F, col.names = T, row.names = F, sep = "\t")

# selecting genomes for CONTROL cohort
#  genome_build %in% GRCh38 , affection_status %in% unaffected, repeat_size < 40, and population %in% EUR since the majority of genomes in the cases cohort correspond to Europeans. 
# Unrelated, belonging to different families each of them.
# reported vs genetic checks PASS
table_ehv3 = read.csv("~/Documents/STRs/ANALYSIS/cases_controls/batch_august/EHv322/table_STR_repeat_size_each_row_allele_EHv3.2.2_HTT_simplified.tsv",
                      stringsAsFactors = F,
                      header = T,
                      sep = "\t")
dim(table_ehv3)
# 310228  27

# Retrieve all expanded genomes for HTT (>=40)
list_expanded = table_ehv3 %>% 
  filter(repeat_size >=40) %>%
  select(platekey) %>%
  pull() %>%
  unique()
length(list_expanded)
# 61

table_ehv3 = left_join(table_ehv3,
                       clin_data %>% select(platekey, family_medical_review_qc_state_code),
                       by = "platekey")

df_controls = table_ehv3 %>%
  filter(!platekey %in% list_expanded, genome_build %in% "GRCh38",
         affection_status %in% "Unaffected",superpopu %in% "EUR",
         family_medical_review_qc_state_code %in% "Passed medical review - for interpretation")
dim(df_controls)
# 562208  28

# Random selection of 30 genomes
set.seed(983823)

# Taking unrelated ones (by FamilyId)
l_random_30_family = sample(df_controls$rare_diseases_family_id, 30)
l_random_100_family = sample(df_controls$rare_diseases_family_id, 100)

random_30_control = df_controls %>%
  filter(rare_diseases_family_id %in% l_random_30_family) %>%
  select(rare_diseases_family_id, platekey, genome_build, participant_phenotypic_sex, superpopu) %>%
  unique()
dim(random_30_control)
# 54  5

random_100_control = df_controls %>%
  filter(rare_diseases_family_id %in% l_random_100_family) %>%
  select(rare_diseases_family_id, platekey, genome_build, participant_phenotypic_sex, superpopu) %>%
  unique()
dim(random_100_control)
# 181  5

# Enrich them with the path to the gVCF
upload_report = read.csv("~/Documents/STRs/data/research/input/batch_august2020_EHv255_and_EHv322/input/upload_report.2020-08-18.txt",
                         stringsAsFactors = F,
                         sep = "\t",
                         comment.char = "#",
                         header = F)
# Remove `_copied` from V3
upload_report$V3 = gsub("_copied", "", upload_report$V3)

random_30_control = left_join(random_30_control,
                              upload_report %>% select(V3,V6),
                              by = c("platekey" = "V3"))
random_100_control = left_join(random_100_control,
                              upload_report %>% select(V3,V6),
                              by = c("platekey" = "V3"))

# Let's finish absolute path of the genome VCF 
random_30_control$V6 = gsub("_copied", "/Variations/", random_30_control$V6)
random_100_control$V6 = gsub("_copied", "/Variations/", random_100_control$V6)

random_30_control = random_30_control %>%
  group_by(platekey) %>%
  mutate(V6 = paste(V6, paste(platekey, ".genome.vcf.gz", sep = ""), sep = "")) %>%
  ungroup() %>%
  as.data.frame()
random_30_control = unique(random_30_control)
dim(random_30_control)
# 54 6

random_100_control = random_100_control %>%
  group_by(platekey) %>%
  mutate(V6 = paste(V6, paste(platekey, ".genome.vcf.gz", sep = ""), sep = "")) %>%
  ungroup() %>%
  as.data.frame()
random_100_control = unique(random_100_control)
dim(random_100_control)
# 195  6

# Write them into files
write.table(random_30_control,
            "./table_30genomes_CONTROL_unrelated_HTT.tsv",
            quote = F,
            row.names = F,
            col.names = T,
            sep = "\t")

# We only need 1 genome from each family
# Creating a new column `elegido` which is a random value from different platekeys from the same family ID
random_100_control = random_100_control %>%
  group_by(rare_diseases_family_id) %>%
  mutate(elegido = sample(platekey, 1)) %>%
  ungroup() %>%
  as.data.frame()

l_selected = unique(random_100_control$elegido)

write.table(random_100_control,
            "./table_100genomes_CONTROL_unrelated_HTT_with_families.tsv",
            quote = F,
            row.names = F,
            col.names = T,
            sep = "\t")

selected_individuals = random_100_control %>%
  filter(platekey %in% l_selected) 

# In selected individuals select the ones that have not been included as Benchmarking
selected_individuals = selected_individuals %>%
  filter(!grepl("BE", V6))

write.table(selected_individuals,
            "./table_99genomes_CONTROL_unrelated_HTT.tsv",
            quote = F,
            row.names = F,
            col.names = T,
            sep = "\t")

# Also write here, the phenotype-part for CONTROLs for later
# ID, IID, sex, CaseControl
pheno_controls = selected_individuals %>% 
  select(platekey, participant_phenotypic_sex)
pheno_controls$ID = seq(from = length(pheno_cases$ID) + 1, to = length(pheno_cases$ID) + length(pheno_controls$platekey), 1)
pheno_controls$ID = paste("FAM", pheno_controls$ID, sep = '_')
pheno_controls$CaseControl = rep("1", length(pheno_controls$platekey))

# format females as `2` and males as `1`
pheno_controls$participant_phenotypic_sex = gsub("Female", "2", pheno_controls$participant_phenotypic_sex)
pheno_controls$participant_phenotypic_sex = gsub("Male", "1", pheno_controls$participant_phenotypic_sex)

# Reorder columns and rename them
pheno_controls = pheno_controls %>%
  select(ID, platekey, participant_phenotypic_sex, CaseControl)
colnames(pheno_controls) = c("ID", "IID", "sex", "CaseControl")
write.table(pheno_controls, "list_phenotypes_controls.tsv", quote = F, col.names = T, row.names = F, sep = "\t")


