# Objective: plot percentage (not frequency) of each repeat-size in 
# - individuals recruited under `neurology`
# - rest of individuals
# - NIID
# - inclusions
# last both sent by Zhongbo
date ()
Sys.info ()[c("nodename", "user")]
commandArgs ()
rm (list = ls ())
R.version.string ## "R version 3.6.2 (2019-12-12)"

# libraries
library(dplyr)
library(ggplot2)

# Set working directory
working_dir="~/Documents/STRs/PAPERS/NIID/"
setwd(working_dir)

# Load data
df_zhongbo = read.csv("data/zhonbo_data_kris_fixed.tsv",
                      sep = "\t",
                      header = T,
                      stringsAsFactors = F)
dim(df_zhongbo)
# 60  2

df_zhongbo$repeat_size = as.integer(df_zhongbo$repeat_size)

# Load data corresponding to 56K genomes (to align with other figures) we do have for populations
popu_table = read.csv("~/Documents/STRs/ANALYSIS/population_research/population_and_super-population_definitions_across_59352_WGS_REv9_271119.tsv",
                               header = T,
                               sep = "\t",
                               stringsAsFactors = F)
dim(popu_table)
# 59356  21

l_genomes = unique(popu_table$platekey)
length(l_genomes)
# 59356

# Retrieve directly from EHv3 merged table the percentages for NOTCH2 gene
notch2_table = read.csv("~/Documents/STRs/ANALYSIS/cases_controls/EHv3.1.2/table_STR_repeat_size_each_row_allele_EHv3.1.2_NOTCH2NL_simplified.tsv",
                      sep = "\t",
                      stringsAsFactors = F,
                      header = T)
dim(notch2_table)
# 152594  19

notch2_table_neuro = notch2_table %>%
  filter(grepl("[Nn][Ee][Uu][Rr][Oo]", disease_group)) 
dim(notch2_table_neuro)
# 24826  19

notch2_table_not_neuro = notch2_table %>%
  filter(!grepl("[Nn][Ee][Uu][Rr][Oo]", disease_group))
dim(notch2_table_not_neuro)
# 127768  19


# Let's create percentage of alleles
notch2_table_neuro = notch2_table_neuro %>%
  group_by(repeat_size) %>%
  mutate(percent_allele = n()/length(notch2_table_neuro$platekey)) %>%
  ungroup() %>%
  as.data.frame()

notch2_table_not_neuro = notch2_table_not_neuro %>%
  group_by(repeat_size) %>%
  mutate(percent_allele = n()/length(notch2_table_not_neuro$platekey)) %>%
  ungroup() %>%
  as.data.frame()


# Adapt zhongbo table in order to have both alleles in the same column
df_zhongbo_NIDD = df_zhongbo %>%
  filter(group %in% "NIID")


df_zhongbo_inclusion = df_zhongbo %>%
  filter(!group %in% "NIID")


df_zhongbo_NIDD = df_zhongbo_NIDD %>%
  group_by(repeat_size) %>%
  mutate(percent_allele = n()/length(df_zhongbo_NIDD$group)) %>%
  ungroup() %>%
  as.data.frame()

df_zhongbo_inclusion = df_zhongbo_inclusion %>%
  group_by(repeat_size) %>%
  mutate(percent_allele = n()/length(df_zhongbo_inclusion$group)) %>%
  ungroup() %>%
  as.data.frame()

notch2_table_neuro = notch2_table_neuro %>%
  select(repeat_size, percent_allele)
notch2_table_neuro$group = rep("Neuro", length(notch2_table_neuro$repeat_size))

notch2_table_not_neuro = notch2_table_not_neuro %>%
  select(repeat_size, percent_allele)
notch2_table_not_neuro$group = rep("Not neuro", length(notch2_table_not_neuro$repeat_size))

merged_all = rbind(df_zhongbo_inclusion,
                   df_zhongbo_NIDD,
                   notch2_table_neuro,
                   notch2_table_not_neuro)

joint_plot = ggplot(merged_all, aes(x = repeat_size, y = percent_allele, group = group, color = group)) + 
  #geom_line() + 
  geom_bar(stat = "identity", position = "dodge") +
  theme_classic() +
  ylab("Allele frequency") + 
  xlab("Number of CGG repeat expansions") 

png("./figures/Figure1A_zhongbo_dodged.png",units="in", width=5, height=5, res=300)
print(joint_plot)
dev.off()


