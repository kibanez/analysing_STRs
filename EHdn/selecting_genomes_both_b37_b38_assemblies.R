# Objective: there are genomes that might have been sequenced in both genome assemblies: b37 and b38
date ()
Sys.info ()[c("nodename", "user")]
commandArgs ()
rm (list = ls ())
R.version.string ## "R version 3.6.1 (2019-07-05)"

# libraries
library(dplyr)

# defining working directory
setwd("~/Documents/STRs/ANALYSIS/EHdn/EHdn-v0.8.6/input/")

# load main clinical data table
main_clin_data = read.csv("~/Documents/STRs/clinical_data/clinical_data/rd_genomes_all_data_041219.tsv",
                          sep = "\t",
                          stringsAsFactors = F,
                          header = T)
dim(main_clin_data)
# 1124633  28


# In the main dataset, we can have genomes sequenced in both genome assemblies: GRCh37 and GRCh38
l_platekey_b37 = main_clin_data %>%
  filter(genome_build %in% "GRCh37") %>%
  select(plate_key) %>%
  unique() %>%
  pull()
length(l_platekey_b37)
# 9901

l_platekey_b38 = main_clin_data %>%
  filter(genome_build %in% "GRCh38") %>%
  select(plate_key) %>%
  unique() %>%
  pull()
length(l_platekey_b38)
# 63551

l_platekey_both = unique(intersect(l_platekey_b37,
                            l_platekey_b38))
length(l_platekey_both)
# 1664

# Taking BAM file in GRCh37 genome assembly
# Let's take the path for all these genomes
# Why? They have been analysed through EHdn with GRCh37 and GRCh38, but the last one was GRCh38...we need to be sure we are using GRCh37 genomes
# We will create 2 files, a list with GRCh37 BAM files (path to the BAM files) and a list with GRCh38 BAM files (path to the BAM files)

df_platekey_both = main_clin_data %>%
  filter(plate_key %in% l_platekey_both, genome_build %in% "GRCh37", grepl("HX", file_path)) %>%
  select(plate_key, file_path)

df_platekey_both = unique(df_platekey_both)
dim(df_platekey_both)
# 1638  2

# There are some that are not captured with `HX`
which_out = setdiff(l_platekey_both, unique(df_platekey_both$plate_key))
df_which_out =  main_clin_data %>%
  filter(plate_key %in% which_out, genome_build %in% "GRCh37") %>%
  select(plate_key, file_path)
df_which_out = unique(df_which_out)

selected = c()
for (i in 1:length(df_which_out$plate_key)){
  l_loci = df_which_out %>% filter(plate_key %in% df_which_out$plate_key[i])
  l_year = unlist(lapply(strsplit(l_loci$file_path, '/'), '[', 4))
  selected = rbind(selected, l_loci[pmatch(max(l_year),l_year),])
}
dim(selected)
# 110  2

selected = unique(selected)
dim(selected)
# 54  2

# Include this to df_platekey_both
df_platekey_both = rbind(df_platekey_both,
                         selected)

# There are still some duplicated genomes with the same genome build (GRCh37)
index_duplicated = which(duplicated(df_platekey_both$plate_key))

# 1- remove them from `df_platekey_both`
l_duplicated_genomes = df_platekey_both$plate_key[index_duplicated]
l_duplicated_genomes = unique(l_duplicated_genomes)
length(l_duplicated_genomes)
# 26

df_platekey_both = df_platekey_both %>%
  filter(!plate_key %in% l_duplicated_genomes)

# 2- select the most recent one
df_platekey_both2 = main_clin_data %>%
  filter(plate_key %in% l_duplicated_genomes) %>%
  select(plate_key, file_path)
df_platekey_both2 = unique(df_platekey_both2)
dim(df_platekey_both2)
# 86  2

selected = c()
for (i in 1:length(df_platekey_both2$plate_key)){
  l_loci = df_platekey_both2 %>% filter(plate_key %in% df_platekey_both2$plate_key[i])
  l_year = unlist(lapply(strsplit(l_loci$file_path, '/'), '[', 4))
  selected = rbind(selected, l_loci[pmatch(max(l_year),l_year),])
}
dim(selected)
# 86  2
selected= unique(selected)
dim(selected)
# 26  2

df_platekey_both = rbind(df_platekey_both,
                         selected)
dim(df_platekey_both)
# 1664  2

df_platekey_both = unique(df_platekey_both)
dim(df_platekey_both)
# 1664  2

length(df_platekey_both$plate_key)
# 1664
length(df_platekey_both$file_path)
# 1664

write.table(df_platekey_both,
            "~/Documents/STRs/ANALYSIS/EHdn/EHdn-v0.8.6/input/list_1664_genomes_both_b37_and_b38_in_GRCh37.csv",
            sep = ",",
            row.names = F,
            col.names = F,
            quote = F)
