# Objective: visualise expanded genomes from batch march 2020, across 13 loci, EHv322
# https://docs.google.com/spreadsheets/d/1cuh2rsDkQP3YEHjX6ogLWlgxO3sd0Jfs4zCVdizBQEc/edit#gid=621033475
date ()
Sys.info ()[c("nodename", "user")]
commandArgs ()
rm (list = ls ())
R.version.string ## "R version 3.6.1 (2019-07-05)"

library(ggplot2); packageDescription ("ggplot2", fields = "Version") #"3.2.1"
library(grid); packageDescription ("grid", fields = "Version") #"3.6.1"
library(gridExtra); packageDescription ("gridExtra", fields = "Version") #"2.3"
library(reshape2); packageDescription ("reshape2", fields = "Version") #"1.4.3"
require(dplyr); packageDescription ("dplyr", fields = "Version") #"0.8.3"

# Set working dir
setwd("~/Documents/STRs/ANALYSIS/population_research/PAPER/expanded_genomes_main_pilot/")

# Load data - all merged
merged_table = read.csv("./batch_march_92K_EHv322_expansions_beyond_premutation - merged_all.tsv",
                        stringsAsFactors = F, 
                        header = T,
                        sep = "\t")
dim(merged_table)
# 514  4

# Load MAIN and PILOT ancestry info, to retrieve PC1 and PC2 values, for each genome to be plotted
popu_table = read.csv("~/Documents/STRs/ANALYSIS/population_research/MAIN_ANCESTRY/GEL_60k_germline_dataset_fine_grained_population_assignment20200224.csv",
                      stringsAsFactors = F, 
                      sep = ",",
                      header = T)
dim(popu_table)
# 59464  36

pilot_popu_table = read.csv("~/Documents/STRs/ANALYSIS/population_research/PILOT_ANCESTRY/FINE_GRAINED_RF_classifications_incl_superPOP_prediction_final20191216.csv",
                            stringsAsFactors = F,
                            sep = ",",
                            header = T)
dim(pilot_popu_table)
# 4821  44 

# recode PC names in pilot_popu_table
colnames(pilot_popu_table) = c("ID", "PC1", "PC2", "PC3", "PC4", "PC5", "PC6", colnames(pilot_popu_table)[c(8:44)])

# We cannot Merge Main and Pilot ancestry tables, since PCs have been computed in a diff way

# enrich merged_table with PC1 and PC2 values
list_exp_genomes = unique(merged_table$ID)
length(list_exp_genomes)
# 512

merged_table_main = left_join(merged_table,
                              popu_table %>% filter(ID %in% list_exp_genomes) %>% select(ID, PC1, PC2),
                              by = "ID")
dim(merged_table_main)
# 514  6

merged_table_pilot = left_join(merged_table,
                              pilot_popu_table %>% filter(ID %in% list_exp_genomes) %>% select(ID, PC1, PC2),
                              by = "ID")
dim(merged_table_pilot)
# 514  6

png("figures/expanded_genomes_MAIN.png")
ggplot(data=merged_table_main %>% filter(!is.na(merged.superpopu)), 
       aes(x=PC2, y=PC1, colour = merged.superpopu)) +
  geom_point() +
  xlab("PC2") +
  ylab("PC1") +
  guides(fill = FALSE)
dev.off()

# only unrelated
merged_table_main_unrelated = merged_table_main %>%
  filter(!is.na(merged.familyID), !is.na(PC1))

l_dup_families = merged_table_main_unrelated$merged.familyID[which(duplicated(merged_table_main_unrelated$merged.familyID))]
# Take only 1 member for each of the duplicated families

dup_families = merged_table_main_unrelated %>%
  filter(merged.familyID %in% l_dup_families)

merged_table_main_unrelated = merged_table_main_unrelated %>%
  filter(!merged.familyID %in% l_dup_families)
dim(merged_table_main_unrelated)
# 200 6

to_add= data.frame()
for(i in 1:length(l_dup_families)){
  aux = dup_families %>%
    filter(merged.familyID %in% l_dup_families[i])
  
  to_add = rbind(to_add,
                 aux[1,])
}

merged_table_main_unrelated = rbind(merged_table_main_unrelated,
                                    to_add)


png("figures/expanded_unrelated_genomes_MAIN.png")
ggplot(data=merged_table_main_unrelated %>% filter(!is.na(merged.superpopu)), 
       aes(x=PC2, y=PC1, colour = merged.superpopu)) +
  geom_point() +
  xlab("PC2") +
  ylab("PC1") +
  guides(fill = FALSE)
dev.off()



# PILOT

png("figures/expanded_unrelated_genomes_PILOT.png")
ggplot(data=merged_table_pilot %>% filter(!is.na(merged.superpopu), !is.na(PC1)), 
       aes(x=PC2, y=PC1, colour = merged.superpopu)) +
  geom_point() +
  xlab("PC2") +
  ylab("PC1") +
  guides(fill = FALSE)
dev.off()


# only unrelated
merged_table_pilot_unrelated = merged_table_pilot %>%
  filter(!is.na(merged.familyID), !is.na(PC1))
dim(merged_table_pilot_unrelated)
# 45  6

l_dup_families = merged_table_pilot_unrelated$merged.familyID[which(duplicated(merged_table_pilot_unrelated$merged.familyID))]
length(l_dup_families)
# 1

# Take only 1 member for each of the duplicated families
dup_families = merged_table_pilot_unrelated %>%
  filter(merged.familyID %in% l_dup_families)

#dup_families
#locus                ID merged.superpopu merged.familyID         PC1         PC2
#1   HTT LP2000711-DNA_B12              EUR        50001165 0.002179575 -0.02190111
#2   FXN LP2000711-DNA_B12              EUR        50001165 0.002179575 -0.02190111
# Diff loci, we are good then



# Barplots

# Let's plot the raw numbers of each ancestry superpopu
raw_numbers_popus_main = as.data.frame(table(merged_table_main_unrelated$merged.superpopu))
colnames(raw_numbers_popus_main) = c("population", "Number of genomes")


raw_numbers_popus_pilot = as.data.frame(table(merged_table_pilot_unrelated$merged.superpopu))
colnames(raw_numbers_popus_pilot) = c("population", "Number of genomes")


png("figures/barplot_ancestry_groups_raw_numbers_main_unrelated.png")
ggplot(raw_numbers_popus_main, 
       aes(x = reorder(population, -`Number of genomes`), y = `Number of genomes`)) + 
  geom_bar(stat = "identity", aes(fill = population)) + 
  geom_text(aes(label=`Number of genomes`), vjust=-0.5, size = 4, colour = "grey") +
  ylab("Number of genomes") + 
  xlab("Ancestry cohorts - EHv3.2.2 - Main cohort") 
dev.off()

png("figures/barplot_ancestry_groups_raw_numbers_pilot_unrelated.png")
ggplot(raw_numbers_popus_pilot, 
       aes(x = reorder(population, -`Number of genomes`), y = `Number of genomes`)) + 
  geom_bar(stat = "identity", aes(fill = population)) + 
  geom_text(aes(label=`Number of genomes`), vjust=-0.5, size = 4, colour = "grey") +
  ylab("Number of genomes") + 
  xlab("Ancestry cohorts - EHv3.2.2 - Pilot cohort") 
dev.off()

raw_numbers_popus_merged = full_join(raw_numbers_popus_main,
                                     raw_numbers_popus_pilot,
                                     by = "population")
# replace NA by 0
raw_numbers_popus_merged$`Number of genomes.y`[which(is.na(raw_numbers_popus_merged$`Number of genomes.y`))] = 0

raw_numbers_popus_merged = raw_numbers_popus_merged %>%
  mutate(number_genomes = `Number of genomes.x` + `Number of genomes.y`)


png("figures/barplot_ancestry_groups_raw_numbers_pilotANDmain_unrelated.png")
ggplot(raw_numbers_popus_merged, 
       aes(x = reorder(population, -number_genomes), y = number_genomes)) + 
  geom_bar(stat = "identity", aes(fill = population)) + 
  geom_text(aes(label=number_genomes), vjust=-0.5, size = 4, colour = "grey") +
  ylab("Number of genomes") + 
  xlab("Ancestry cohorts - EHv3.2.2") 
dev.off()

# Who are the 15 AFR genomes?
# Load data - all merged
merged_table = read.csv("./batch_march_92K_EHv322_expansions_beyond_premutation - merged_all.tsv",
                        stringsAsFactors = F, 
                        header = T,
                        sep = "\t")
dim(merged_table)
# 514  4

to_print_AFR = merged_table %>% filter(merged.superpopu %in% "AFR", !is.na(merged.familyID))
write.table(to_print_AFR,
            "table_17_related_genomes_AFR.tsv",
            quote = F,
            row.names = F,
            col.names = T,
            sep = "\t")





