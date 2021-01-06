# Objective: analyse haplotyping blocks within HTT
date ()
Sys.info ()[c("nodename", "user")]
commandArgs ()
rm (list = ls ())
R.version.string ## "R version 3.6.1 (2019-07-05)"

# libraries
library(dplyr)

# Set working directory
setwd("~/Documents/STRs/ANALYSIS/haplotyping/HTT/")

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

# Retrieve genme build, sex, popu for these genomes
haplo_genomes = clin_data %>%
  filter(platekey %in% l_exp_genomes) %>%
  select(platekey, genome_build, participant_ethnic_category, participant_phenotypic_sex, superpopu, programme, file_path)
haplo_genomes = unique(haplo_genomes)
dim(haplo_genomes)
# 29  7

length(unique(haplo_genomes$platekey))
# 29

# write them into a file
write.table(haplo_genomes,
            "./table_29genomes_expanded_unrelated_HTT.tsv",
            sep = "\t",
            quote = F,
            row.names = F,
            col.names = T)

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


