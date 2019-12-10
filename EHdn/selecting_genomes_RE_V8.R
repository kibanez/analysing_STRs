# Objective: generate list of genomes to run EHdn through: germline RD and cancer genomes

date ()
Sys.info ()[c("nodename", "user")]
commandArgs ()
rm (list = ls ())
R.version.string ## "R version 3.6.1 (2019-07-05)"

# libraries
library(dplyr)
library(ggplot2)
library(reshape)

setwd("/Users/kibanez/Documents/STRs/ANALYSIS/EHdn/")

# Let's take the germline RD and cancer table from RE V8
germline_table = read.csv("~/Documents/STRs/clinical_data/clinical_data/rd_genomes_all_data_041219.tsv",
                          sep = "\t",
                          header = T,
                          stringsAsFactors = F)
dim(germline_table)
# 1124633  28

all_data_V8 = read.csv("~/Documents/STRs/clinical_data/clinical_data/raw/RE_clinical_data_V8/genome_file_paths_and_types_2019-12-04_15-13-29.tsv",
                sep = "\t",
                stringsAsFactors = F,
                header=T)
dim(all_data_V8)
#cancer germline        cancer somatic experimental germline  experimental somatic rare disease germline 
#73813                 52249                   268                   159                373954 

# Let's take germlineL cancer and RD

subset_V8 = all_data_V8 %>%
  filter(grepl("germline", type), file_sub_type %in% "BAM") %>%
  select(platekey, file_path, genome_build)
dim(subset_V8)
# 90252  3

length(unique(subset_V8$platekey))
# 88086

# Write into a file
write.table(subset_V8, 
            "./input_EHdn-v0.8.6/list_90252_genomes_88086_unique_RE_V8_all.csv",
            sep = ",",
            quote = F,
            row.names = F,
            col.names = F)

# Subset of GRCh37
write.table(subset_V8 %>% filter(genome_build %in% "GRCh37") %>% select(platekey, file_path), 
            "input_EHdn-v0.8.6/list_90252_genomes_88086_unique_RE_V8_b37.csv",
            sep = ",",
            quote = F,
            row.names = F,
            col.names = F)

# Subset of GRCh38
write.table(subset_V8 %>% filter(genome_build %in% "GRCh38") %>% select(platekey, file_path), 
            "input_EHdn-v0.8.6/list_90252_genomes_88086_unique_RE_V8_b38.csv",
            sep = ",",
            quote = F,
            row.names = F,
            col.names = F)
