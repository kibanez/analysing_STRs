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
setwd("/Users/kibanez/Documents/STRs/ANALYSIS/population_research/PAPER/expanded_genomes_main_pilot/unrelated/")

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

popu_table_unrel = popu_table %>%
  filter(ID %in% l_unrelated)
dim(popu_table_unrel)
# 38344 36

# Enrich with superpopu
popu_table_unrel = popu_table_unrel %>%
  mutate(superpopu = case_when(best_guess_predicted_ancstry == "PJL" ~ "SAS",
                               best_guess_predicted_ancstry == "GBR" ~ "EUR",
                               best_guess_predicted_ancstry == "CEU" ~ "EUR",
                               best_guess_predicted_ancstry == "TSI" ~ "EUR",
                               best_guess_predicted_ancstry == "PUR" ~ "AMR",
                               best_guess_predicted_ancstry == "ACB" ~ "AFR",
                               best_guess_predicted_ancstry == "GIH" ~ "SAS",
                               best_guess_predicted_ancstry == "ASW" ~ "AFR",
                               best_guess_predicted_ancstry == "MXL" ~ "AMR",
                               best_guess_predicted_ancstry == "ESN" ~ "AFR",
                               best_guess_predicted_ancstry == "LWK" ~ "AFR",
                               best_guess_predicted_ancstry == "CHS" ~ "EAS",
                               best_guess_predicted_ancstry == "BEB" ~ "SAS",
                               best_guess_predicted_ancstry == "KHV" ~ "EAS",
                               best_guess_predicted_ancstry == "CLM" ~ "AMR",
                               best_guess_predicted_ancstry == "MSL" ~ "AFR",
                               best_guess_predicted_ancstry == "YRI" ~ "AFR",
                               best_guess_predicted_ancstry == "GWD" ~ "AFR",
                               best_guess_predicted_ancstry == "FIN" ~ "EUR",
                               best_guess_predicted_ancstry == "ITU" ~ "SAS",
                               best_guess_predicted_ancstry == "JPT" ~ "EAS",
                               best_guess_predicted_ancstry == "STU" ~ "SAS",
                               best_guess_predicted_ancstry == "CHB" ~ "EAS",
                               best_guess_predicted_ancstry == "PEL" ~ "AMR",
                               best_guess_predicted_ancstry == "IBS" ~ "EUR"))

# Load here individual tables, per locus, including expanded genomes -- we need to focus only in unrelated genomes!
l_files = list.files("./PAPER/expanded_genomes_main_pilot/", pattern = "\\.tsv$")
l_files = paste("./PAPER/expanded_genomes_main_pilot/", l_files, sep = "")

l_loci = c("AR", "ATN1", "ATXN1", "ATXN2", "ATXN3", "ATXN7", "CACNA1A", "C9ORF72", "DMPK", "HTT", "FMR1", "FXN", "TBP")
l_superpopu = c("AFR", "AMR", "EAS", "EUR", "SAS")
for(i in 1:length(l_files)){
  print(l_files[i])
  locus_table = read.csv(l_files[i],
                         stringsAsFactors = F,
                         header = T,
                         sep = "\t")
  
  # focus on unrelated genomes
  locus_table = locus_table %>%
    filter(ID %in% l_unrelated)
  
  # Print total number of unrelated genomes with an expansion
  locus_table %>% select(ID) %>% unique() %>% pull() %>% length() %>% print()
  
  # How many genomes in each superpopu?
  for (j in 1:length(l_superpopu)){
    print(l_superpopu[j])
    
    num_superpopu = locus_table %>%
      filter(merged_superpopu %in% l_superpopu[j]) %>%
      select(ID) %>%
      unique() %>%
      pull() %>%
      length() %>%
      print()
    
  }
}


atn1_table_unrel = atn1_table %>% filter(ID %in% l_unrelated) 
dim(atn1_table_unrel)
# 11 27


# With all population coctails
png("figures/population_distribution_unrelated_all_ancestries.png")
ggplot(data=popu_table_unrel %>% filter(!is.na(self_reported)), 
       aes(x=PC2, y=PC1, colour = superpopu)) +
  geom_point() +
  xlab("PC2 across 38,344 unrelated genomes") +
  ylab("PC1 across 38,344 unrelated genomes") +
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