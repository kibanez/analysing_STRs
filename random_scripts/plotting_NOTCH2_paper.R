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

# For Supp Figure2, we want to focus on NEURO patients: families that have been recruited under `NEURO` as `disease_group`
l_family_neuro = notch2_table %>%
  filter(grepl("[Nn][Ee][Uu][Rr][Oo]", disease_group)) %>%
  select(rare_diseases_family_id) %>%
  unique() %>%
  pull() 
length(l_family_neuro)
# 10920

# Retrieve all genomes within this list
notch2_table_neuro = notch2_table %>%
  filter(rare_diseases_family_id %in% l_family_neuro)
dim(notch2_table_neuro)
# 49876  19

# Now, we only can plot those genomes for which we do have ancestry info
l_participant_notch2_neuro = unique(notch2_table_neuro$participant_id)
length(l_participant_notch2_neuro)
# 24206

length(unique(notch2_table_neuro$platekey))
# 24313

# how many genomes/platekeys for only super ancestries??
notch2_table_neuro %>% 
  filter(population %in% c("AFR", "AMR", "ASI", "EAS", "EUR")) %>%
  select(platekey) %>%
  unique() %>%
  pull() %>%
  length()
# 19699


# There is no need to use `popu_table` since I already did this for case-control tables !! :)
violin_popu_plot = ggplot(notch2_table_neuro %>% filter(population %in% c("AFR", "AMR", "ASI", "EAS", "EUR")), 
                          aes(x = population, y = repeat_size, fill = population)) + 
  geom_violin() +
  geom_boxplot(width=0.1, fill="white") +
  theme_classic() +
  ylab("number of repeats") + 
  xlab("ethnicity") 

png("./figures/SuppFigure2.png",units="in", width=5, height=5, res=600)
print(violin_popu_plot)
dev.off()


# Table with IQR
# Focusing on those having super-population info (the ones in the violin plot)
notch2_table_neuro_superpopu = notch2_table_neuro %>% 
  filter(population %in% c("AFR", "AMR", "ASI", "EAS", "EUR")) 

# Compute how many platekeys for each super population and IQR for repeat-sizes
# EUR
notch2_table_neuro_superpopu %>% filter(population %in% "EUR") %>% select(platekey) %>% unique() %>% pull() %>% length()
# 16443
notch2_table_neuro_superpopu %>% filter(population %in% "EUR") %>% select(repeat_size) %>% summary()
#repeat_size  
#Min.   : 5.0  
#1st Qu.:15.0  
#Median :20.0  
#Mean   :19.2  
#3rd Qu.:22.0  
#Max.   :60.0  

# EAS
notch2_table_neuro_superpopu %>% filter(population %in% "EAS") %>% select(platekey) %>% unique() %>% pull() %>% length()
# 110
notch2_table_neuro_superpopu %>% filter(population %in% "EAS") %>% select(repeat_size) %>% summary()
#repeat_size   
#Min.   :10.00  
#1st Qu.:15.00  
#Median :20.00  
#Mean   :20.67  
#3rd Qu.:24.00  
#Max.   :55.00


# AMR
notch2_table_neuro_superpopu %>% filter(population %in% "AMR") %>% select(platekey) %>% unique() %>% pull() %>% length()
# 315
notch2_table_neuro_superpopu %>% filter(population %in% "AMR") %>% select(repeat_size) %>% summary()
#repeat_size   
#Min.   : 7.00  
#1st Qu.:16.75  
#Median :20.00  
#Mean   :19.55  
#3rd Qu.:22.00  
#Max.   :56.00 


# Figure 1 (not anymore)
# Let's create percentage of alleles
notch2_table_neuro = notch2_table_neuro %>%
  group_by(repeat_size) %>%
  mutate(percent_allele = n()/length(notch2_table_neuro$platekey)) %>%
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

png("./figures/Figure1A_zhongbo_dodged.png", units="in", width=5, height=5, res=300)
print(joint_plot)
dev.off()


