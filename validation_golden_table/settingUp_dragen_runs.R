# Objective: setup the input file with the CRAM files we do have frmo dragen v3.2
# All have been realigned to GRCh38, in CRAM format
# The format we want to have is
# <PLATEKEY> \t <PATH> \t <GENDER> \t <GENOME_BUILD>
date ()
Sys.info ()[c("nodename", "user")]
commandArgs ()
rm (list = ls ())
R.version.string ## "R version 3.6.1 (2019-07-05)"

# libraries
library(dplyr)
library(tidyr)

# working directory
setwd("~/Documents/STRs/VALIDATION/dragen/input")

all_data = read.csv("list_fam_platekey_path_dragenv3.2.csv",
                    sep = ",",
                    stringsAsFactors = F,
                    header = F)

clin_data = read.csv("~/Documents/STRs/clinical_data/clinical_data/rd_genomes_all_data_041219.tsv",
                     sep = "\t",
                     header = T,
                     stringsAsFactors = F)
pilot_clin_data = read.csv("~/Documents/STRs/clinical_data/pilot_clinical_data/pilot_cohort_clinical_data_4833_genomes_removingPanels_280919.tsv",
                           sep = "\t",
                           stringsAsFactors = F,
                           header = T)

# We need to run EHv2/EHv3 across genomes included in our validation dataset
nhnn_data = read.table("NHNN_negative_families_b37.txt", stringsAsFactors = F)$V1
ncl_data = read.table("NCL_negative_famIDs_b37.txt", stringsAsFactors = F)$V1
wessex_data1 = read.table("Wessex_negative_famIDs.txt", stringsAsFactors = F)$V1
wessex_data2 = read.table("Wessex_negative_families_b37.txt", stringsAsFactors = F)$V1
others= read.table("all_STR_fams_neg_and_pos.txt", stringsAsFactors = F)$V1
others2 = read.table("str_ISO_families.txt", stringsAsFactors = F)$V1

l_genomes = c(nhnn_data,
              ncl_data,
              wessex_data1,
              wessex_data2,
              others,
              others2)

length(l_genomes)
# 220

# Retrieve from all_data the info about these 220 genomes

val_data = all_data %>% filter(V1 %in% l_genomes)
dim(val_data)
# 245  3

colnames(val_data) = c("familyID", "platekey", "path")

# Let's enrich val_data with gender
val_data = left_join(val_data,
                     clin_data %>% select(platekey, participant_phenotypic_sex))

val_data = unique(val_data)
dim(val_data)
# 245  4

val_data = left_join(val_data,
                     pilot_clin_data %>% select(plateKey, sex),
                     by = c("platekey" = "plateKey"))


# There are some pilot...
l_na = which(is.na(val_data$participant_phenotypic_sex))

val_data$participant_phenotypic_sex[l_na] = val_data$sex[l_na]

val_data = unique(val_data)
dim(val_data)
# 245  5

# Write into the final table
write.table(val_data %>% select(platekey, path, participant_phenotypic_sex),
            "list_genomes_to_runEH_from_dragen_alignments.tsv",
            sep = "\t",
            quote = F,
            row.names = F,
            col.names = F)

# Check how many genomes withint the golden validation set are included in this run Dorota did with dragen v3.2
list_unique_platekeys_golden_table = read.table("~/Documents/STRs/VALIDATION/list_platekeys_unique.tsv",
                                                stringsAsFactors = F)
list_unique_platekeys_golden_table = list_unique_platekeys_golden_table$V1

# Intersection
length(intersect(list_unique_platekeys_golden_table, unique(val_data$platekey)))
# 93

# Which are the ones not included there
l_genomes_extra_to_run_dragen = setdiff(list_unique_platekeys_golden_table,
                                        unique(val_data$platekey))



  
  




