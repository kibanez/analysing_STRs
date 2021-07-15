# Objective: given a list of platekeys, we want to take the repeat-sizes for genes
date ()
Sys.info ()[c("nodename", "user")]
commandArgs ()
rm (list = ls ())
R.version.string ## "R version 3.3.2 (2016-10-31)"

# libraries
library(dplyr); packageDescription ("dplyr", fields = "Version") #"0.8.5"
library(tidyverse); packageDescription ("tidyverse", fields = "Version") # "1.2.1
library(ggplot2); packageDescription ("ggplot2", fields = "Version") #"3.3.0"

# set the working directory
setwd("~/Downloads/")

# Load the data
df_mito = read.csv("~/Downloads/table_mito_genomes.tsv",
                   stringsAsFactors = F,
                   header = T,
                   sep = "\t")
dim(df_mito)
# 345 2

# load merged august batch 2020
merged_table = read.csv("~/Documents/STRs/data/research/batch_august2020/output_EHv3.2.2_vcfs/merged/merged_93446_genomes_EHv322_batch_august2020.tsv",
                        stringsAsFactors = F,
                        header = T,
                        sep = "\t")
dim(merged_table)
# 27238  12

# Can you see if they are in your data for the population paper, and pull out a the repeat-sizes for all the 13 RE loci for all genomes? 
# A table with rows= genomes and columns =  26 (13 loci x 2 alleles)
list_genes = c("AR", "ATN1", "ATXN1", "ATXN2", "ATXN3", "ATXN7", "C9ORF72", "CACNA1A", "DMPK", "FMR1", "FXN","HTT", "TBP")
merged_table = merged_table %>%
  filter(gene %in% list_genes)

# Load list of unrelated genomes
l_unrel = read.csv("~/Documents/STRs/ANALYSIS/population_research/MAIN_ANCESTRY/batch2/l_unrelated_55603_genomes_batch2.txt",
                   stringsAsFactors = F,
                   header = F)
l_unrel = l_unrel$V1
length(l_unrel)
# 55603

l_platekeys = df_mito$platekey
df_genomes = data.frame()
for(i in 1:length(l_platekeys)){
  is_unrel = l_platekeys[i] %in% l_unrel
  merged_platekey = merged_table %>%
    filter(grepl(l_platekeys[i], list_samples))
  
   for (j in 1:length(list_genes)){
     merged_platekey_locus = merged_platekey %>%
       filter(gene %in% list_genes[j]) 
     if (dim(merged_platekey_locus)[1] > 1){
       to_extract = merged_platekey_locus %>% select(gene, allele)
       min_EHv3 = min(to_extract$allele)
       max_EHv3 = max(to_extract$allele)
       
       df_genomes = rbind(df_genomes,
                          cbind(l_platekeys[i], list_genes[j], min_EHv3, max_EHv3, is_unrel))
     }else if(dim(merged_platekey_locus)[1] == 0){
       min_EHv3 = "lowQ"
       max_EHv3 = "lowQ"
       df_genomes = rbind(df_genomes,
                          cbind(l_platekeys[i], list_genes[j], min_EHv3, max_EHv3, is_unrel))
     }else{
       # check whether only has 1 allele (XY sample and AR or FMR1 gene) or it's duplicated (_x2)
       if (grepl(paste(l_platekeys[i], ".vcf_x2", sep = ""), merged_platekey_locus$list_samples)){
         to_extract = merged_platekey_locus %>% select(gene, allele)
         min_EHv3 = to_extract$allele
         max_EHv3 = to_extract$allele
         
         df_genomes = rbind(df_genomes,
                            cbind(l_platekeys[i], list_genes[j], min_EHv3, max_EHv3, is_unrel))
       }else{
         to_extract = merged_platekey_locus %>% select(gene, allele)
         min_EHv3 = to_extract$allele
         max_EHv3 = "NA"
         
         df_genomes = rbind(df_genomes,
                            cbind(l_platekeys[i], list_genes[j], min_EHv3, max_EHv3, is_unrel))
       }
     }
   }
}
dim(df_genomes)
# 4485  5

colnames(df_genomes) = c("platekey", "gene", "min_EHv3", "max_EHv3", "is_unrel")
write.table(df_genomes, "~/Downloads/table_mito_genomes.tsv", quote = F, col.names = F, row.names = F, sep = "\t")