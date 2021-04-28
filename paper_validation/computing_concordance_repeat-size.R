# Objective: compute repeat-size concordance between EH and PCR, just keeping 1 gene for those duplicated or repeated genomes
date ()
Sys.info ()[c("nodename", "user")]
commandArgs ()
rm (list = ls ())
R.version.string ## "R version 3.3.2 (2016-10-31)"

# libraries
library(dplyr); packageDescription ("dplyr", fields = "Version") #"0.8.5"
library(ggplot2); packageDescription ("ggplot2", fields = "Version") #"3.3.0"

# set the working directory
setwd("~/Documents/STRs/PAPERS/VALIDATION_PAPER/")

concordance = read.csv("./concordance_repeats.tsv", stringsAsFactors = F, sep = "\t", header = T)
dim(concordance)
# 902  11

platekey_duplicated = unique(concordance$platekey[which(duplicated(concordance$platekey))])
length(platekey_duplicated)
# 195

# There are a total of 195 repeat platekeys - for which we've got PCR lengths for more than one gene
# We first take all dedup platekeys
concordance_dedup = concordance %>%
  filter(!platekey %in% platekey_duplicated) %>%
  select(locus, validation_id, platekey, PCR, EHv312, concordance, classi, type)

# We take only 1 gene/locus for each repeated platekey 
l_genes = unique(concordance$locus)
set.seed(532)

concordance_dup = concordance %>%
  filter(platekey %in% platekey_duplicated) %>%
  group_by(platekey) %>%
  mutate(random_gene = sample(l_genes, 1)) %>%
  ungroup() %>%
  as.data.frame()

concordance_dup = concordance_dup %>%
  select(random_gene, validation_id, platekey, PCR, EHv312, concordance, classi, type)
colnames(concordance_dup) = colnames(concordance_dedup)
concordance_dup = unique(concordance_dup)
dim(concordance_dup)
# 737  8

concordance_merge = rbind(concordance_dedup,
                          concordance_dup)
write.table(concordance_merge,
          "./concordance_repeats_dedup_repeated_platekeys.tsv",
          quote = F, row.names = F, col.names = T, sep = "\t")
