# Objective: analyse the repeat-size distribution across 9 genes shared by Sharp et al, after running EHdn on Parkinsonism
# For EHv322
date ()
Sys.info ()[c("nodename", "user")]
commandArgs ()
rm (list = ls ())
R.version.string ## "R version 3.6.1 (2019-07-05)"

# libraries
library(ggplot2); packageDescription ("ggplot2", fields = "Version") #"3.2.1"
library(grid); packageDescription ("grid", fields = "Version") #"3.6.1"
library(gridExtra); packageDescription ("gridExtra", fields = "Version") #"2.3"
library(reshape2); packageDescription ("reshape2", fields = "Version") #"1.4.3"
require(dplyr); packageDescription ("dplyr", fields = "Version") #"0.8.3"
library(magick); packageDescription ("magick", fields = "Version") #"2.4.0"
library(cowplot); packageDescription ("cowplot", fields = "Version") #"1.0.0"
library(ggpubr); packageDescription ("ggpubr", fields = "Version") #"1.0.0"
library(purrr); packageDescription ("purrr", fields = "Version") #"0.3.3"

# Set working dir
setwd("~/Documents/STRs/ANALYSIS/SHARP/EHdn_Parkinson/analysis/")

# Functions
source("~/git/analysing_STRs/functions/plot_together_histo_boxplot.R")
source("~/git/analysing_STRs/functions/computing_percentiles.R")

# load merged august data
merged_data = read.csv("~/Documents/STRs/ANALYSIS/SHARP/EHdn_Parkinson/output_EHv3.2.2_vcfs/merged/merged_93425_genomes_EHv322.tsv",
                       stringsAsFactors = F,
                       sep = "\t",
                       header = T)
dim(merged_data)
# 795  12

# load clinical data - changing to RE V11 (since we are sharing with external groups)
clin_data = read.csv("~/Documents/STRs/clinical_data/clinical_data/Main_RE_V11_and_Pilot_programmes.tsv",
                     sep = "\t",
                     stringsAsFactors = F,
                     header = T)
dim(clin_data)
# 2444984  26

# Only focusing on germline genomes
clin_data = clin_data %>%
  filter(grepl("germline", type))

# Load type (case/control/pseudocase/pseudocontrol) for genomes
type_data = read.csv("./table_cases_controls_84518_genomes_cases_controls_pseudoca_pseudoco.csv",
                     stringsAsFactors = F)
dim(type_data)
# 84518  2

clin_data = left_join(clin_data,
                      type_data,
                      by = "platekey")
dim(clin_data)
# 2417423  27

# List of cases/controls/pseudocases/pseudocontrols
l_cases = clin_data %>%
  filter(type.y %in% "case") %>%
  select(platekey) %>%
  unique() %>%
  pull()
length(l_cases)
# 14478

l_controls = clin_data %>%
  filter(type.y %in% "control") %>%
  select(platekey) %>%
  unique() %>%
  pull()
length(l_controls)
# 24095

l_pseudocases = clin_data %>%
  filter(type.y %in% "pseudocase") %>%
  select(platekey) %>%
  unique() %>%
  pull()
length(l_pseudocases)
# 16857

l_pseudocontrols = clin_data %>%
  filter(type.y %in% "pseudocontrol") %>%
  select(platekey) %>%
  unique() %>%
  pull()
length(l_pseudocontrols)
# 13600

# 1. Merge GRCh37 and GRCh38 info, since chromosome names are different
# GRCh38 are chr1, chr2, chr3 while GRCh37 are 1,2,3
# In SHARP everything should be GRCh38, because Andy sent us coordinates only in GRCh38
merged_data$chr = recode(merged_data$chr,
                         "1" = "chr1",
                         "2" = "chr2",
                         "3" = "chr3",
                         "4" = "chr4",
                         "5" = "chr5",
                         "6" = "chr6",
                         "7" = "chr7",
                         "8" = "chr8",
                         "9" = "chr9",
                         "10" = "chr10",
                         "11" = "chr11",
                         "12" = "chr12",
                         "13" = "chr13",
                         "14" = "chr14",
                         "15" = "chr15",
                         "16" = "chr16",
                         "17" = "chr17",
                         "18" = "chr18",
                         "19" = "chr19",
                         "20" = "chr20",
                         "21" = "chr21",
                         "22" = "chr22",
                         "X" = "chrX")

# Let's focus on SHARP genes
l_genes = unique(merged_data$gene)
length(l_genes)
# 329

# Let's focus only on SHARP genes
l_sharp = l_genes[which(grepl("SHARP", l_genes, ignore.case = TRUE))]
length(l_sharp)
# 55

sharp_merged_data = merged_data %>%
  filter(gene %in% l_sharp)
dim(sharp_merged_data)
# 5974  12

# Output folder
output_folder = 'EHv322_batch_august2020'
#dir.create(output_folder)

for(i in 1:length(l_sharp)){
  plot_together_histo_boxplot(df_input = sharp_merged_data,
                              gene_name = l_sharp[i],
                              output_folder = output_folder,
                              l_platekeys_probands_neuro_unique,
                              l_platekeys_probands_notNeuro_unique)
}

# Summarise report with quantiles for all genes in SHARP
# Based on 93k genomes
df_percentiles = computing_percentiles(sharp_merged_data)

# Write down first a line with total number of probands in neuro and probands in not neuro
to_write = cbind("93446",
                 length(l_platekeys_probands_neuro_unique),
                 length(l_platekeys_probands_notNeuro_unique))
to_write = rbind(c("Total number of genomes","Total number of probands recruited under Neurological disorders", "Total number of probands NOT recruited under Neurological disorders"),
                 to_write,
                 cbind("", "", ""))
write.table(to_write, 
            "./EHv322_batch_august2020/summary_stats_quantiles_55_SHARP_genes.tsv", 
            sep='\t', 
            quote = F,
            row.names=F, 
            col.names=F)
write.table(df_percentiles,
            "./EHv322_batch_august2020/summary_stats_quantiles_55_SHARP_genes.tsv", 
            append = T, 
            sep = "\t", 
            quote = F,
            row.names = F,
            col.names = T)
