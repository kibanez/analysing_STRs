# Objective: we want to analyse the populations that are enriched in the pathogenic tails for each locus within 1Kg Phase3 dataset
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
setwd("~/Documents/STRs/ANALYSIS/population_research/PAPER/expanded_genomes_main_pilot/1Kg/")

# load data
merged_table = read.csv("~/Documents/STRs/ANALYSIS/population_research/1kg/data/all_data/merged/merged_all_2504_genomes_1Kg_EHv3.2.2.tsv",
                        sep = "\t",
                        stringsAsFactors = F, 
                        header = T)
dim(merged_table)
# 1342  12

# Load 1Kg metadata
metadata = read.csv("~/Documents/STRs/ANALYSIS/population_research/1kg/1000G_2504_high_coverage.sequence.index.tsv",
                    stringsAsFactors = F,
                    sep = "\t",
                    header = T)
dim(metadata)
# 2504  22

metadata = metadata %>% select(SAMPLE_NAME, POPULATION)
dim(metadata)
# 2504  2

# Recode superpopu
metadata = metadata %>%
  mutate(superpopu = case_when(POPULATION == "PJL" ~ "SAS",
                               POPULATION == "GBR" ~ "EUR",
                               POPULATION == "CEU" ~ "EUR",
                               POPULATION == "TSI" ~ "EUR",
                               POPULATION == "PUR" ~ "AMR",
                               POPULATION == "ACB" ~ "AFR",
                               POPULATION == "GIH" ~ "SAS",
                               POPULATION == "ASW" ~ "AFR",
                               POPULATION == "MXL" ~ "AMR",
                               POPULATION == "ESN" ~ "AFR",
                               POPULATION == "LWK" ~ "AFR",
                               POPULATION == "CHS" ~ "EAS",
                               POPULATION == "BEB" ~ "SAS",
                               POPULATION == "KHV" ~ "EAS",
                               POPULATION == "CLM" ~ "AMR",
                               POPULATION == "MSL" ~ "AFR",
                               POPULATION == "YRI" ~ "AFR",
                               POPULATION == "GWD" ~ "AFR",
                               POPULATION == "FIN" ~ "EUR",
                               POPULATION == "ITU" ~ "SAS",
                               POPULATION == "JPT" ~ "EAS",
                               POPULATION == "STU" ~ "SAS",
                               POPULATION == "CHB" ~ "EAS",
                               POPULATION == "CDX" ~ "EAS",
                               POPULATION == "PEL" ~ "AMR",
                               POPULATION == "IBS" ~ "EUR"))

dim(metadata)
# 2504  3

# For each locus
l_genes = c("AR", "ATN1", "ATXN1", "ATXN2", "ATXN3", "ATXN7", "CACNA1A", "C9ORF72", "DMPK", "HTT", "FMR1", "FXN", "TBP")
l_premut_cutoff = c(34,34,35,31,43,17,17,30,50,35,55,44,41)

for (i in 1:length(l_genes)){
  print(l_genes[i])
  merged_table_locus = merged_table %>%
    filter(gene %in% l_genes[i], allele > l_premut_cutoff[i])
  
  if (dim(merged_table_locus)[1] > 0){
    list_allele_size = c()
    list_vcf_patho_locus = c()
    df_platekey_size = data.frame()
    for (j in 1:length(merged_table_locus$list_samples)){
      list_vcf_allele = strsplit(merged_table_locus$list_samples[j], ';')[[1]]
      number_vcf = length(list_vcf_allele)
      list_vcf_patho_locus = c(list_vcf_patho_locus,
                               list_vcf_allele)
      list_allele_size = rep(merged_table_locus$allele[j], number_vcf)
      
      list_vcf_allele = gsub('.vcf', '', list_vcf_allele)
      list_vcf_allele = gsub('^EHv3.2.2_', '', list_vcf_allele)
      list_vcf_allele = gsub('_x2', '', list_vcf_allele)
      
      # Create dataframe with platekey-repeat-size, for the expanded genomes
      df_platekey_size = rbind(df_platekey_size,
                               data.frame(platekey = list_vcf_allele,
                                          repeat_size = list_allele_size))
      df_platekey_size$platekey = as.character(df_platekey_size$platekey)
      df_platekey_size$repeat_size = as.integer(df_platekey_size$repeat_size)
    }
    
    list_vcf_patho_locus = gsub('.vcf', '', list_vcf_patho_locus)
    list_vcf_patho_locus = gsub('^EHv3.2.2_', '', list_vcf_patho_locus)
    list_vcf_patho_locus = gsub('_x2', '', list_vcf_patho_locus)
    
    # Enrich platekeys now with ancestry info
    patho_popu = metadata %>%
      filter(SAMPLE_NAME %in% list_vcf_patho_locus) %>%
      select(SAMPLE_NAME, POPULATION, superpopu)
    print(dim(patho_popu))
    
    # Add locus name as column
    patho_popu$locus = rep(l_genes[i], length(patho_popu$SAMPLE_NAME))
    
    # merge
    patho_merged = left_join(df_platekey_size,
                             patho_popu,
                             by = c("platekey" = "SAMPLE_NAME"))
    
    output_file_name = paste(l_genes[i], "beyond_", sep = "_")
    output_file_name = paste(output_file_name, "premutation_cutoff_", sep = "_")
    output_file_name = paste(output_file_name, as.character(l_premut_cutoff[i]), sep = "")
    output_file_name = paste(output_file_name, "EHv322_1Kg.tsv", sep = "_")
    write.table(patho_merged, 
                output_file_name, 
                sep = "\t",
                quote = F,
                row.names = F,
                col.names = T)
  }
}


# we have merged in google excel all loci
merged_expanded = read.csv("./1K_phase3_EHv322_expansions_beyond_premutation - merged.tsv",
                           stringsAsFactors = F,
                           header = T,
                           sep = "\t")
dim(merged_expanded)  
# 109  5

