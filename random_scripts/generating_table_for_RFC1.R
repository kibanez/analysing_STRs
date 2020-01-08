# Objective: generate table with <PLATEKEY> | <VCF_FILE> | <PATH> | <BUILD>
date ()
Sys.info ()[c("nodename", "user")]
commandArgs ()
rm (list = ls ())
R.version.string ## "R version 3.6.1 (2019-07-05)"

# libraries
library(dplyr)

# Set working directory
setwd("~/Documents/STRs/ANALYSIS/CANVAS_RFC1/")

# Upload latest clinical data
clin_data = read.csv("~/Documents/STRs/clinical_data/clinical_data/rd_genomes_all_data_041219.tsv",
                     stringsAsFactors = F,
                     sep = "\t",
                     header = TRUE)
dim(clin_data)
# 1124633  28

clin_data = clin_data %>% select(platekey, file_path, genome_build)
clin_data = unique(clin_data)
dim(clin_data)
# 111148  3


# Upload the list of genomes for EHv2.5.5
list_ehv2 = read.table("list_7159_genomes_RFC1_EHv2.5.5.tsv", stringsAsFactors = F)
list_ehv2 = list_ehv2$V1
length(list_ehv2)
# 7159

list_ehv2 = gsub('.{4}$', '', list_ehv2)
list_ehv2 = gsub('^.{3}', '', list_ehv2)

clin_data_ehv2 = clin_data %>% 
  filter(platekey %in% list_ehv2)
dim(clin_data_ehv2)
# 7548  2

clin_data_ehv2 = unique(clin_data_ehv2)
dim(clin_data_ehv2)

write.table(clin_data_ehv2,
            "table_7548_genomes_for_7159_unique_genomes_RFC1_EHv2.5.5.tsv",
            quote = F,
            row.names = F,
            col.names = T,
            sep = "\t")