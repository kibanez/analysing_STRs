# Objective: from EHv2.5.5 summer batch select genomes with the following selection criteria (and compare with Arianna's numbers):
# both from GRCH37 and GRCH38:
# proband only
# specific disease == intellectual disability
# Arianna's numbers are the following:
# 5783 probands recruited under ID
# repeat size > 55
# 131 probands
date ()
Sys.info ()[c("nodename", "user")]
commandArgs ()
rm (list = ls ())
R.version.string ## "R version 3.6.1 (2019-07-05)"

library(ggplot2); packageDescription ("ggplot2", fields = "Version") #"3.2.1"
library(grid); packageDescription ("grid", fields = "Version") #"3.6.2"
library(gridExtra); packageDescription ("gridExtra", fields = "Version") #"2.3"
library(reshape2); packageDescription ("reshape2", fields = "Version") #"1.4.3"
require(dplyr); packageDescription ("dplyr", fields = "Version") #"0.8.3"

# Set working dir
setwd("~/Documents/STRs/data/research/EH_2.5.5_research_August2019/EH_output_v2.5.5_research_August_2019/FMR1_GRCh37_GRCh38_encrypted/")

# Load clinically enriched data for FMR1
table_b37 = read.csv("./research_genomes_86457_GRCh37_EHv2.5.5_FMR1.tsv",
                        sep = "\t",
                        stringsAsFactors = F,
                        header = T)
dim(table_b37)
# 11913  26

table_b38 = read.csv("./research_genomes_86457_GRCh38_EHv2.5.5_FMR1_CGG.tsv",
                     sep = "\t",
                     stringsAsFactors = F,
                     header = T)
dim(table_b38)
# 108939  26

# merge b37 and b38 tables
merged_table = rbind(table_b37,
                     table_b38)

merged_table = unique(merged_table)
dim(merged_table)
# 120852  26

# Recode b37 and b38 chrX name
merged_table$chr = recode(merged_table$chr,
                          "X" = "chrX")

# 1- take only probands (== "N/A) [is.na() are genomes from Cancer programme]
fmr1_filtering = merged_table %>%
  filter(biological_relationship_to_proband %in% "N/A")
length(unique(fmr1_filtering$participant_id))
# 26530

# 2 - `specific disease` == ID
fmr1_filtering = fmr1_filtering %>%
  filter(grepl("[Ii]ntellectual disability",specific_disease))
length(unique(fmr1_filtering$participant_id))
# 5037

# 3 - repeat.size > 55
fmr1_filtering = fmr1_filtering %>%
  filter(repeat.size > 55)
length(unique(fmr1_filtering$participant_id))
# 130

# Arianna's work is ending with 131 probands. 
# Let's see which 2 are different
ari_fmr1 = read.table("~/Documents/STRs/VALIDATION/for_ILM_pileup/list_131_FMR1_probands_Arianna.txt", stringsAsFactors = F)
ari_fmr1 = ari_fmr1$V1
length(ari_fmr1)
# 131

l_130_my_analysis = unique(fmr1_filtering$plate_key.x)
length(l_130_my_analysis)
# 130

# differences?
setdiff(ari_fmr1, l_130_my_analysis)
# "LP3000046-DNA_E01" 

# let's see what they are

merged_table %>%
  filter(plate_key.x %in% "LP3000046-DNA_E01") %>%
  View()

  



