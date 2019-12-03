# Objective: creating tables for each locus, in order to play with cases/controls and see possible existing signals
# For EHv2.5.5
date ()
Sys.info ()[c("nodename", "user")]
commandArgs ()
rm (list = ls ())
R.version.string ## "R version 3.6.1 (2019-07-05)"

# libraries
library(dplyr)

# Load data
merged_data = read.csv("/Users/kibanez/Documents/STRs/data/research/EH_2.5.5_research_August2019/EH_output_v2.5.5_research_August_2019/merged_genotypeUpdated/merged_loci_86457_research_genomes_new_loci_EHv2.5.5_summer2019_removingListVcfFiles.tsv",
                       sep = '\t',
                       header = T,
                       stringsAsFactors = F)
dim(merged_data)
# 3983  11

# As output
# We want to the following output - NOTE each row is an allele (!!!)
# Family ID
# Patient ID
# LP number
# Participant type (proband / mother / brother etc)
# HTT allele
# Disease Category
# Disease Subgroup
# Panels applied
# Sex
# Age at recruitment
# Age of onset
# no of family participants
# ethnic 

l_genes = unique(merged_data$gene)

#for (){
 i = "HTT_CAG" 
 locus_data = merged_data %>% filter(gene %in% i)
 
#}
