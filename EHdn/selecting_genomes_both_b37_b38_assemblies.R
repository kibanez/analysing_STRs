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

write.table(l_platekey_both,
            "~/Documents/STRs/ANALYSIS/EHdn/EHdn-v0.8.6/input/list_1664_genomes_both_b37_and_b38.txt",
            quote = F,
            row.names = F,
            col.names = F)

df_platekey_both = main_clin_data %>%
  filter(plate_key %in% l_platekey_both, genome_build %in% "GRCh37", grepl("HX", file_path)) %>%
  select(plate_key, file_path)

df_platekey_both = unique(df_platekey_both)
dim(df_platekey_both)
# 1638  2

# There are still some duplicated genomes with the same genome build (GRCh37)
index_duplicated = which(duplicated(df_platekey_both$plate_key))
#[1]   92   94   95  114  118  327  534  548  950 1013 1016 1018 1036 1043 1052 1059 1067 1115 1294 1307 1309 1311 1312 1314 1402
#[26] 1410 1588 1604 1681 1682

# 1- remove them from `df_platekey_both`
l_duplicated_genomes = df_platekey_both$plate_key[index_duplicated]
df_platekey_both = df_platekey_both[-index_duplicated,]

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
  selected = rbind(selected, l_loci[which(max(l_year) %in% l_year),])
}

df_platekey_both = rbind(df_platekey_both,
                         selected)
dim(df_platekey_both)
# 1696  2



write.table(df_platekey_both,
            "~/Documents/STRs/ANALYSIS/EHdn/EHdn-v0.8.6/input/list_1664_genomes_both_b37_and_b38.csv",
            sep = ",",
            row.names = F,
            col.names = F,
            quote = F)
