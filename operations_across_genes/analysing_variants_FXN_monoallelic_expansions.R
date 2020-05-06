# Objective: there is a list of genomes that have a pathogenic expansion on FXN, but only in one allele. 
# we want to analyse whether they also have a potentially pathogenic variant in FXN
date ()
Sys.info ()[c("nodename", "user")]
commandArgs ()
rm (list = ls ())
R.version.string ## "R version 3.6.3 (2020-02-29)"

# libraries
library(dplyr)

# defining working directory
setwd("~/Documents/STRs/ANALYSIS/FXN_monoallelic/")

# 1 - Retrieve the genome.vcf.gz for each genome
clin_data = read.csv("~/Documents/STRs/clinical_data/clinical_data/rd_genomes_all_data_300320.tsv",
                     sep = "\t",
                     stringsAsFactors = F,
                     header = T)
dim(clin_data)
# 1124633      31

list_fxn_genomes = read.table("./list_43_platekeys_monoallelic_STR_FXN.tsv", stringsAsFactors = F)
list_fxn_genomes = list_fxn_genomes$V1
length(list_fxn_genomes)
# 43

path_vcf = clin_data %>%
  filter(platekey %in% list_fxn_genomes) %>%
  select(file_path) %>%
  unique() %>%
  pull()
length(path_vcf)
# 43

