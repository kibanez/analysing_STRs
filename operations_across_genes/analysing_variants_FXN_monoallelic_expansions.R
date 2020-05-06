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

l_path_genome = c()
for (i in 1:length(path_vcf)){
  l_split = strsplit(path_vcf[i], "/")[[1]]
  # take the part of the path we are interested in
  l_split = l_split[c(1:6)]
  l_split = c(l_split, "Variations")
  l_split = c(l_split, paste(l_split[6], ".genome.vcf.gz", sep = ""))
  new_char = paste(l_split, collapse = '/')
  l_path_genome = c(l_path_genome,
                    new_char)
  
}
length(l_path_genome)
# 43
