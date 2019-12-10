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

# I saw that genomes considered for ancestry computation, were not in my list for germline - cancer in particular
# I want to rescue here genomes in V7 were taken into account for ancestry and we haven't include them in the final list
l_genomes_for_ehdn = unique(subset_V8$platekey)
length(l_genomes_for_ehdn)
# 88086

popu_table = read.csv("~/Documents/STRs/ANALYSIS/population_research/population_and_super-population_definitions_across_59352_WGS_REv9_271119.tsv",
                      header = T,
                      stringsAsFactors = F,
                      sep = "\t")
dim(popu_table)
# 59356  21

length(intersect(l_genomes_for_ehdn, unique(popu_table$platekey)))
# 58441

l_genomes_in_popu_not_ehdn = setdiff(unique(popu_table$platekey), l_genomes_for_ehdn)
length(l_genomes_in_popu_not_ehdn)
# 915

#Let's take their path to the BAM files (from V7 since the popu table is froml_geno  V7)
# ALso, see if there are some germline genomes...not in V8
all_data_V7 = read.csv("~/Documents/STRs/clinical_data/clinical_data/raw/RE_clinical_data_V7/genome_file_paths_and_types_2019-12-05_21-35-16_V7.tsv",
                       sep = "\t",
                       stringsAsFactors = F,
                       header=T)
dim(all_data_V7)
# 492236  10

subset_V7 = all_data_V7 %>%
  filter(grepl("germline", type), file_sub_type %in% "BAM") %>%
  select(platekey, file_path, genome_build)
dim(subset_V7)
# 87323  3

l_genomes_V7 = unique(subset_V7$platekey)
length(l_genomes_V7)
# 85532

length(intersect(l_genomes_V7, l_genomes_for_ehdn))
# 84475

length(setdiff(l_genomes_V7, l_genomes_for_ehdn))
# 1057

length(intersect(l_genomes_V7, l_genomes_in_popu_not_ehdn))
# 562

l_genomes_in_popu_not_ehdn_from_V7 = setdiff(l_genomes_V7, l_genomes_for_ehdn)

l_genomes_new_for_eh = unique(c(l_genomes_in_popu_not_ehdn_from_V7, l_genomes_in_popu_not_ehdn_from_V7))
length(l_genomes_new_for_eh)
# 1057

# Let's take the paths

df_genomes_new_for_eh = subset_V7 %>% 
  filter(platekey %in% l_genomes_new_for_eh) %>%
  select(platekey, file_path, genome_build)

dim(df_genomes_new_for_eh)
# 1116  3

# mmmm some of them are "unknown"...let's see what Bertha says
table(df_genomes_new_for_eh$genome_build)
#GRCh37  GRCh38 unknown 
#353     746      17 
       
# They all are from 2018-02-17 --> GRCh38 then
df_genomes_new_for_eh = df_genomes_new_for_eh %>% 
  mutate(genome_build = replace(genome_build, genome_build=="unknown", "GRCh38")) %>%
  as.data.frame()
dim(df_genomes_new_for_eh)
# 1116  3

# Check again
table(df_genomes_new_for_eh$genome_build)
# GRCh37 GRCh38 
# 353    763

# Write results into table
write.table(df_genomes_new_for_eh,
            "input_EHdn-v0.8.6/list_1057_genomes_rescued_from_RE_V7_all.csv",
            sep = ",",
            quote = F,
            col.names = F,
            row.names = F)

# GRCh37
write.table(df_genomes_new_for_eh %>% filter(genome_build %in% "GRCh37") %>% select(platekey, file_path),
            "input_EHdn-v0.8.6/list_1057_genomes_rescued_from_RE_V7_b37.csv",
            sep = ",",
            quote = F,
            col.names = F,
            row.names = F)


# GRCh38
write.table(df_genomes_new_for_eh %>% filter(genome_build %in% "GRCh38") %>% select(platekey, file_path),
            "input_EHdn-v0.8.6/list_1057_genomes_rescued_from_RE_V7_b38.csv",
            sep = ",",
            quote = F,
            col.names = F,
            row.names = F)






