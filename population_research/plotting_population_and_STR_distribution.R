# Objective: take advantage from the work GEL bioinfo research group has done estimating ancestry for ~59,356 genomes
# and analyse and study the STR distribution across them
date ()
Sys.info ()[c("nodename", "user")]
commandArgs ()
rm (list = ls ())
R.version.string ## "R version 3.6.1 (2019-07-05)"

# libraries
library(dplyr); packageDescription ("dplyr", fields = "Version") #"0.8.3"
library(ggplot2); packageDescription ("ggplot2", fields = "Version") #"3.2.1"
library(reshape); packageDescription ("reshape", fields = "Version") #"0.8.8"
library(scatterplot3d); packageDescription("scatterplot3d", fields = "Version") # 0.3-41
library(ggpubr); packageDescription("ggpubr", fields = "Version") # 0.2.3


# Set environment
setwd("/Users/kibanez/Documents/STRs/ANALYSIS/population_research/")

# Load population data
popu_table = read.csv("~/Documents/STRs/ANALYSIS/population_research/MAIN_ANCESTRY/GEL_60k_germline_dataset_fine_grained_population_assignment20200224.csv",
                      stringsAsFactors = F,
                      header = T)
dim(popu_table)
# 59464  36

popu_table_pilot = read.csv("~/Documents/STRs/ANALYSIS/population_research/PILOT_ANCESTRY/FINE_GRAINED_RF_classifications_incl_superPOP_prediction_final20191216.csv",
                            stringsAsFactors = F,
                            header = T)
dim(popu_table_pilot)
# 4821  44

# Load thresholds
# STR annotation, threshold including the largest normal and the smallest pathogenic sizes reported
gene_annotation_normal = '/Users/kibanez/git/analysing_STRs/threshold_largest_normal_reported_research.txt'
gene_data_normal = read.table(gene_annotation_normal, stringsAsFactors=F, header = T)

gene_annotation_pathogenic = '/Users/kibanez/git/analysing_STRs/threshold_smallest_pathogenic_reported_research.txt'
gene_data_pathogenic = read.table(gene_annotation_pathogenic, stringsAsFactors=F, header = T)


# Functions
source("/Users/kibanez/git/analysing_STRs/functions/plot_violin_ancestry.R")
source("/Users/kibanez/git/analysing_STRs/functions/plot_gene.R")
source("/Users/kibanez/git/analysing_STRs/functions/plot_gene_joint_ancestries.R")

# Let's first explore the data and see the population distribution across 59K
popu_table_enriched = popu_table_enriched %>%
  mutate(population = case_when(pred_african_ancestries >= 0.7 ~ "AFR",
                                pred_american_ancestries >= 0.7 ~ "AMR",
                                pred_european_ancestries >= 0.7 ~ "EUR",
                                pred_east_asian_ancestries >= 0.7 ~ "EAS",
                                pred_south_asian_ancestries >= 0.7 ~ "ASI",
                                pred_african_ancestries >= 0.25 & pred_american_ancestries >= 0.25  ~ "AFR-AMR",
                                pred_african_ancestries >= 0.25 & pred_european_ancestries >= 0.25 ~ "AFR-EUR",
                                pred_african_ancestries >= 0.25 & pred_east_asian_ancestries >= 0.25  ~ "AFR-EAS",
                                pred_african_ancestries >= 0.25 & pred_south_asian_ancestries >= 0.25  ~ "AFR-ASI",
                                pred_american_ancestries >= 0.25 & pred_european_ancestries >= 0.25  ~ "AMR-EUR",
                                pred_american_ancestries >= 0.25 & pred_east_asian_ancestries >= 0.25  ~ "AMR-EAS",
                                pred_american_ancestries >= 0.25 & pred_south_asian_ancestries >= 0.25  ~ "AMR-ASI",
                                pred_european_ancestries >= 0.25 & pred_american_ancestries >= 0.25  ~ "EUR-AMR",
                                pred_european_ancestries >= 0.25 & pred_east_asian_ancestries >= 0.25  ~ "EUR-EAS",
                                pred_european_ancestries >= 0.25 & pred_south_asian_ancestries >= 0.25  ~ "EUR-ASI",
                                pred_east_asian_ancestries >= 0.25 & pred_south_asian_ancestries >= 0.25 ~ "EAS-ASI"))

# Let's write this info into a file, in order to have super-population and population info / genome
write.table(popu_table_enriched,
            "./population_and_super-population_definitions_across_59352_WGS_REv9_271119.tsv",
            quote = F,
            row.names = F,
            col.names = T,
            sep = "\t")

# With all population coctails
png("figures/population_distribution_56176_all_ancestries.png")
ggplot(data=popu_table_enriched %>% filter(!is.na(population)), 
       aes(x=pc2, y=pc1, colour = population)) +
  geom_hex(bins=300) +
  xlab("PC2 across 56,176 genomes") +
  ylab("PC1 across 56,176 genomes") +
  guides(fill = FALSE)
dev.off()

# Just "pure" ancestries
popu_table_enriched %>% filter(population %in% c("AFR", "EUR", "AMR", "EAS", "ASI")) %>% select(platekey) %>% unique() %>% pull() %>% length()
# 56176

png("figures/population_distribution_56176_pure_ancestries.png")
ggplot(data=popu_table_enriched %>% filter(population %in% c("AFR", "EUR", "AMR", "EAS", "ASI")) , 
       aes(x=pc2, y=pc1, colour = population)) +
  geom_hex(bins=300) +
  xlab("PC2 across 56,176 genomes") +
  ylab("PC1 across 56,176 genomes") +
  guides(fill = FALSE)
dev.off()

# Let's plot the raw numbers of each ancestry sub-cohort or sub-group
raw_numbers_popus = as.data.frame(table(popu_table_enriched$population))
colnames(raw_numbers_popus) = c("population", "Number of genomes")

raw_numbers_popus_pure = raw_numbers_popus %>% 
  filter(population %in% c("AFR", "EUR", "AMR", "EAS", "ASI")) 

# Pure ancestries
png("figures/barplot_pure_ancestry_groups_raw_numbers.png")
ggplot(raw_numbers_popus_pure, 
       aes(x = reorder(population, -`Number of genomes`), y = `Number of genomes`)) + 
  geom_bar(stat = "identity", aes(fill = population)) + 
  geom_text(aes(label=`Number of genomes`), vjust=-0.5, size = 4, colour = "grey") +
  ylab("Number of genomes") + 
  xlab("Ancestry cohorts - EHv3.1.2 - 56,176 total genomes") 
dev.off()

# All ancestries
png("figures/barplot_all_ancestry_groups_raw_numbers.png")
ggplot(raw_numbers_popus, 
       aes(x = reorder(population, -`Number of genomes`), y = `Number of genomes`)) + 
  geom_bar(stat = "identity", aes(fill = population)) + 
  geom_text(aes(label=`Number of genomes`), vjust=-0.5, size = 4, colour = "grey") +
  ylab("Number of genomes") + 
  xlab("Ancestry cohorts - EHv3.1.2 - 56,176 total genomes") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
dev.off()


# After plotting the ancestry distribution, let's see how different locus-STR-distribution change across different population cohorts
# let's take only the genomes for which we do have popu information 
l_genomes = popu_table_enriched %>%
  filter(population %in% c("AFR", "EUR", "AMR", "EAS", "ASI")) %>%
  select(platekey) %>%
  unique() %>%
  pull()

length(l_genomes)
# 56176

# Write into a file the list of platekeys for which we do have population info, in order to create the corresponding merge VCF files
write.table(l_genomes,
            file = "list_genomes_56176_pure_ancestry_info.txt",
            quote = F,
            row.names = F,
            col.names = F)

# Let's stratify genomes by ancestry
l_genomes_AFR = popu_table_enriched %>%
  filter(population %in% 'AFR') %>%
  select(platekey) %>%
  unique() %>%
  pull()
length(l_genomes_AFR)
# 1794

l_genomes_EUR = popu_table_enriched %>%
  filter(population %in% 'EUR') %>%
  select(platekey) %>%
  unique() %>%
  pull()
length(l_genomes_EUR)
# 47208

l_genomes_AMR = popu_table_enriched %>%
  filter(population %in% 'AMR') %>%
  select(platekey) %>%
  unique() %>%
  pull()
length(l_genomes_AMR)
# 803

l_genomes_EAS = popu_table_enriched %>%
  filter(population %in% 'EAS') %>%
  select(platekey) %>%
  unique() %>%
  pull()
length(l_genomes_EAS)
# 400

l_genomes_ASI = popu_table_enriched %>%
  filter(population %in% 'ASI') %>%
  select(platekey) %>%
  unique() %>%
  pull()
length(l_genomes_ASI)
# 5971

write.table(l_genomes_AFR,
            file = "list_genomes_1794_pure_AFR.txt",
            quote = F,
            row.names = F,
            col.names = F)
write.table(l_genomes_EUR,
            file = "list_genomes_47208_pure_EUR.txt",
            quote = F,
            row.names = F,
            col.names = F)
write.table(l_genomes_AMR,
            file = "list_genomes_803_pure_AMR.txt",
            quote = F,
            row.names = F,
            col.names = F)
write.table(l_genomes_EAS,
            file = "list_genomes_400_pure_EAS.txt",
            quote = F,
            row.names = F,
            col.names = F)
write.table(l_genomes_ASI,
            file = "list_genomes_5971_pure_ASI.txt",
            quote = F,
            row.names = F,
            col.names = F)


# We then copy these VCF files, and merge with the python script
# Load EH STR output data
# Research ~80K genomes, EH-v3.1.2- November 2019 -- but the ones that are ALSO included in the population info table
df_asi = read.csv('~/Documents/STRs/ANALYSIS/population_research/EH_3.1.2_research_October2019_55419_genomes/ASI/merged/merged_population_genomes_5789_avg_EHv3.1.2_ASI.tsv',
              sep = '\t',
              stringsAsFactors = F,
              header = T)
dim(df_asi)
# 1648 12

df_eas = read.csv('~/Documents/STRs/ANALYSIS/population_research/EH_3.1.2_research_October2019_55419_genomes/EAS/merged/merged_population_genomes_392_avg_EHv3.1.2_EAS.tsv',
                  sep = '\t',
                  stringsAsFactors = F,
                  header = T) 
dim(df_eas)
# 806  12

df_afr = read.csv('~/Documents/STRs/ANALYSIS/population_research/EH_3.1.2_research_October2019_55419_genomes/AFR/merged/merged_population_genomes_1743_avg_EHv3.1.2_AFR.tsv',
                  sep = '\t',
                  stringsAsFactors = F,
                  header = T) 
dim(df_afr)
# 1271  12

df_amr = read.csv('~/Documents/STRs/ANALYSIS/population_research/EH_3.1.2_research_October2019_55419_genomes/AMR/merged/merged_population_genomes_788_avg_EHv3.1.2_AMR.tsv',
                  sep = '\t',
                  stringsAsFactors = F,
                  header = T) 
dim(df_amr)
# 1141  12

df_eur = read.csv('~/Documents/STRs/ANALYSIS/population_research/EH_3.1.2_research_October2019_55419_genomes/EUR/merged/merged_population_genomes_45690_avg_EHv3.1.2_EUR.tsv',
                  sep = '\t',
                  stringsAsFactors = F,
                  header = T) 
dim(df_eur)
# 2424  12


df_asi = df_asi %>% mutate(population = "ASI")
df_eas = df_eas %>% mutate(population = "EAS")
df_afr = df_afr %>% mutate(population = "AFR")
df_amr = df_amr %>% mutate(population = "AMR")
df_eur = df_eur %>% mutate(population = "EUR")

df_all = rbind(df_afr,
               df_amr,
               df_eur,
               df_eas,
               df_asi)
dim(df_all)
# 7290  13

# Population enriched genomes are only GRCh38, we will ignore then GRCh37
output_folder = "./figures/"

l_loci = sort(unique(df_all$gene))
for (i in 1:length(l_loci)){
  # Each locus - Individually
  plot_gene(df_afr, l_loci[i], gene_data_normal, gene_data_pathogenic, output_folder, "GRCh38", "AFR")
  plot_gene(df_amr, l_loci[i], gene_data_normal, gene_data_pathogenic, output_folder, "GRCh38", "AMR")
  plot_gene(df_eur, l_loci[i], gene_data_normal, gene_data_pathogenic, output_folder, "GRCh38", "EUR")
  plot_gene(df_eas, l_loci[i], gene_data_normal, gene_data_pathogenic, output_folder, "GRCh38", "EAS")
  plot_gene(df_asi, l_loci[i], gene_data_normal, gene_data_pathogenic, output_folder, "GRCh38", "ASI")
  
  # Jointly - distribution
  plot_gene_joint_ancestries(df_all, l_loci[i], gene_data_normal, gene_data_pathogenic, output_folder)
  
  # Jointly - Violin plots
  plot_violin_ancestry(df_all, l_loci[i], gene_data_normal, gene_data_pathogenic, output_folder)
}