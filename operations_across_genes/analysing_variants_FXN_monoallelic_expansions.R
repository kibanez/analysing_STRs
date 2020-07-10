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

df_path_vcf = clin_data %>%
  filter(platekey %in% list_fxn_genomes) %>%
  select(file_path, genome_build) %>%
  unique()
dim(df_path_vcf)  
# 43  2

table(df_path_vcf$genome_build)
#GRCh37 GRCh38 
#5     38 

# we need to split into GRCh37 and GRCh38 !!!! in order to annotate them properly

path_vcf_38 = clin_data %>%
  filter(platekey %in% list_fxn_genomes, genome_build %in% "GRCh38") %>%
  select(file_path) %>%
  unique() %>%
  pull()
length(path_vcf_38)
# 38

path_vcf_37 = clin_data %>%
  filter(platekey %in% list_fxn_genomes, genome_build %in% "GRCh37") %>%
  select(file_path) %>%
  unique() %>%
  pull()
length(path_vcf_37)
# 5


l_path_genome38 = c()
for (i in 1:length(path_vcf_38)){
  l_split = strsplit(path_vcf_38[i], "/")[[1]]
  # take the part of the path we are interested in
  l_split = l_split[c(1:6)]
  l_split = c(l_split, "Variations")
  l_split = c(l_split, paste(l_split[6], ".genome.vcf.gz", sep = ""))
  new_char = paste(l_split, collapse = '/')
  l_path_genome38 = c(l_path_genome38,
                      new_char)
  
}
length(l_path_genome38)
# 38

l_path_genome37 = c()
for (i in 1:length(path_vcf_37)){
  l_split = strsplit(path_vcf_37[i], "/")[[1]]
  # take the part of the path we are interested in
  l_split = l_split[c(1:6)]
  l_split = c(l_split, "Variations")
  l_split = c(l_split, paste(l_split[6], ".genome.vcf.gz", sep = ""))
  new_char = paste(l_split, collapse = '/')
  l_path_genome37 = c(l_path_genome37,
                      new_char)
  
}
length(l_path_genome37)
# 5

write.table(l_path_genome37, "list_5_genomeVCF_path_GRCh37.tsv", quote = F, row.names = F, col.names = F)
write.table(l_path_genome38, "list_38_genomeVCF_path_GRCh38.tsv", quote = F, row.names = F, col.names = F)


# July 2020

# defining working directory
setwd("~/Documents/STRs/ANALYSIS/FXN_monoallelic/list_july/")
list_fxn_genomes = read.table("./FXN_single_expansion_list", stringsAsFactors = F, header = T)
list_fxn_genomes = list_fxn_genomes$platekey
length(list_fxn_genomes)
# 1220

# Let's load batch march input data
df_genomes_b37 = read.csv("~/Documents/STRs/data/research/input/batch_march2020_EHv2.5.5_and_EHv3.2.2/input/list_13024_ouf_of_92669_genomes_GRCh37.csv",
                          stringsAsFactors = F,
                          header = T,
                          sep = ",")
dim(df_genomes_b37)
# 13023  3

df_genomes_b38 = read.csv("~/Documents/STRs/data/research/input/batch_march2020_EHv2.5.5_and_EHv3.2.2/input/list_79645_ouf_of_92669_genomes_GRCh38.csv",
                          stringsAsFactors = F,
                          header = T,
                          sep = ",")
dim(df_genomes_b38)
# 79644  3

df_path_vcf = clin_data %>%
  filter(platekey %in% list_fxn_genomes) %>%
  select(file_path, genome_build) %>%
  unique()
dim(df_path_vcf)  
# 1183  2



table(df_path_vcf$genome_build)
#GRCh37 GRCh38 
#129     835
