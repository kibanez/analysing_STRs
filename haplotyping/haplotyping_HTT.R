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
clin_data = read.csv("~/Documents/STRs/clinical_data/clinical_data/Main_RE_V10_and_Pilot_programmes.tsv",
                     sep = "\t",
                     stringsAsFactors = F,
                     header = T)
dim(clin_data)
# 2101385  24

# Load unrel genomes
l_unrel = read.table("~/Documents/STRs/ANALYSIS/population_research/MAIN_ANCESTRY/batch2/l_unrelated_55603_genomes_batch2.txt",
                     stringsAsFactors = F)
l_unrel = l_unrel$V1
length(l_unrel)
# 55603

# Load 29 expanded genomes after visual inspection and unrelated
l_exp_genomes = read.table("./list_28_platekeys_expanded_cases.txt", stringsAsFactors = F)
l_exp_genomes = l_exp_genomes$V1
length(l_exp_genomes)
# 28

# Retrieve genome build, sex, popu for these genomes
haplo_genomes = clin_data %>%
  filter(platekey %in% l_exp_genomes) %>%
  select(platekey, genome_build, participant_ethnic_category, participant_phenotypic_sex, superpopu, programme)
haplo_genomes = unique(haplo_genomes)
dim(haplo_genomes)
# 28  6

length(unique(haplo_genomes$platekey))
# 28

# write them into a file
write.table(haplo_genomes,
            "./table_28genomes_expanded_unrelated_HTT.tsv",
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
# Unrelated, batch2 (aggV2)
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

# Only proband is enriched with genetic_vs_reported checks. Let's take the famIDs that have been undergone successfully
# Only rare diseases. We don't do this with cancer
l_famID_passesGvsRChecks = table_ehv3 %>%
  filter(genetic_vs_reported_results %in% "familyPassesGvsRChecks") %>%
  select(rare_diseases_family_id) %>%
  unique() %>%
  pull()
length(l_famID_passesGvsRChecks)
# 29925

df_controls = table_ehv3 %>%
  filter(!platekey %in% list_expanded, genome_build %in% "GRCh38",
         affection_status %in% "Unaffected",superpopu %in% "EUR",
         rare_diseases_family_id %in% l_famID_passesGvsRChecks)
dim(df_controls)
# 72754  27

# Select only unrelated genomes
df_controls = df_controls %>%
  filter(platekey %in% l_unrel)
df_controls = unique(df_controls)
dim(df_controls)
# 56770 27

l_controls_unrel = unique(df_controls$platekey)
length(l_controls_unrel)
# 18739

# Enrich them with the absolute path to gVCF file
upload_report = read.csv("~/Documents/STRs/clinical_data/clinical_data/upload_report.2020-08-18.txt",
                         stringsAsFactors = F,
                         sep = "\t",
                         header = T)
dim(upload_report)
# 120711  10

# we need to remove `_copied` - IT people changed again
upload_report$Platekey = gsub("_copied", "", upload_report$Platekey)
upload_report$Path = gsub("_copied", "", upload_report$Path)

# create new path
upload_report = upload_report %>%
  group_by(Platekey) %>%
  mutate(gvcf_path = paste(paste(paste(paste(Path, "Variations", sep = "/")),Platekey,sep = "/"),".genome.vcf.gz", sep = "")) %>%
  ungroup() %>%
  as.data.frame()

# Some platekeys have been sequenced more than once, let's select the latest one
upload_report = upload_report %>%
  group_by(Platekey) %>%
  mutate(latest_gvcf_path = max(gvcf_path)) %>%
  ungroup() %>%
  as.data.frame()

# Retrive the gVCF files for the cases
list_gvcf_cases = upload_report %>%
  filter(Platekey %in% l_exp_genomes) %>%
  select(latest_gvcf_path) %>%
  unique() %>%
  pull()
length(list_gvcf_cases)
# 26

# Write into file
write.table(list_gvcf_cases, "list_28_gVCF_controls_unrel_GRCh38_EUR_HTT.txt", quote = F, row.names = F, col.names = F)

# Retrieve the gVCF files for all controls
list_gvcf_controls = upload_report %>%
  filter(Platekey %in% l_controls_unrel, Delivery.Version %in% "V4") %>%
  select(latest_gvcf_path) %>%
  unique() %>%
  pull()
length(list_gvcf_controls)
# 17878

# Write into file
write.table(list_gvcf_controls, "list_17878_gVCF_controls_unrel_GRCh38_EUR_HTT.txt", quote = F, row.names = F, col.names = F)

# Create phenotype file for Plink
list_total_samples = read.table("./cases_and_controls/list_28_cases_17878_controls.txt", stringsAsFactors = F)
list_total_samples = list_total_samples$V1
length(list_total_samples)
# 17906

# FID IID PID MID sex affection
df_phenotype = clin_data %>%
  filter(platekey %in% list_total_samples) %>%
  select(platekey, participant_phenotypic_sex) %>%
  unique()
length(unique(df_phenotype$platekey))
# 17906

# Enrich with affection status following plink format
#-9 missing 
#0 missing
#1 unaffected
#2 affected
df_phenotype = df_phenotype %>%
  group_by(platekey) %>%
  mutate(affection = ifelse(platekey %in% l_exp_genomes, "2", "1")) %>%
  ungroup() %>%
  as.data.frame()

# Reformat gender (male to 1, female to 2)
df_phenotype$participant_phenotypic_sex = gsub("Male", "1", df_phenotype$participant_phenotypic_sex)
df_phenotype$participant_phenotypic_sex = gsub("Female", "2", df_phenotype$participant_phenotypic_sex)

# Include FID, PID, MID
df_phenotype$FID = df_phenotype$platekey
df_phenotype$PID = rep(0, length(df_phenotype$platekey))
df_phenotype$MID = rep(0, length(df_phenotype$platekey))




