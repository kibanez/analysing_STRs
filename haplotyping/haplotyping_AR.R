# Objective: from Matteo/Ari work analysing expansions on AR (from EHv255)
date ()
Sys.info ()[c("nodename", "user")]
commandArgs ()
rm (list = ls ())
R.version.string ## "R version 3.6.1 (2019-07-05)"

# libraries
library(dplyr)

# Set working directory
setwd("~/Documents/STRs/ANALYSIS/haplotyping/AR/HaploView/")

# Load clinical data
clin_data = read.csv("~/Documents/STRs/clinical_data/clinical_data/rd_genomes_all_data_251120_V10.tsv",
                     sep = "\t",
                     stringsAsFactors = F,
                     header = T)
dim(clin_data)
# 2096500  36

# Load 38 expanded genomes in AR- after visual QC inspection
l_exp_genomes = read.table("./list_38_genomes_beyond_patho.txt", stringsAsFactors = F)
l_exp_genomes = l_exp_genomes$V1
length(l_exp_genomes)
# 38 

# Retrieve genme build, sex, popu for these genomes
haplo_genomes = clin_data %>%
  filter(platekey %in% l_exp_genomes) %>%
  select(platekey, genome_build, participant_ethnic_category, participant_phenotypic_sex, superpopu)
haplo_genomes = unique(haplo_genomes)
dim(haplo_genomes)
# 41  5 

# There is 1 genome not included in the clinical data table, which belongs to data release V9
setdiff(l_exp_genomes,unique(haplo_genomes$platekey))
# "LP3000458-DNA_D11"

to_add = c("LP3000458-DNA_D11", "GRCh38", "White: British", "Male", "EUR")
haplo_genomes = rbind(haplo_genomes,
                      to_add)
dim(haplo_genomes)
# 42  5 

# Create 2 groups: Females and Males
haplo_genomes_female = haplo_genomes %>%
  filter(participant_phenotypic_sex %in% "Female")
  
haplo_genomes_male = haplo_genomes %>%
  filter(participant_phenotypic_sex %in% "Male")


length(unique(haplo_genomes_female$platekey))
# 18
length(unique(haplo_genomes_male$platekey))
# 20

write.table(haplo_genomes_female,
            "./table_18genomes_female.tsv",
            sep = "\t",
            quote = F,
            row.names = F,
            col.names = T)

write.table(haplo_genomes_male,
            "./table_20genomes_male.tsv",
            sep = "\t",
            quote = F,
            row.names = F,
            col.names = T)
