# Objective: create a single table including all 13 loci, and 90,863 genomes for which we have repeat-size for all these genes
# Enrich with relatedness (`Yes` or `No`)
# Enrich with population (batch1 AND batch2)
# For EHv3.2.2
date ()
Sys.info ()[c("nodename", "user")]
commandArgs ()
rm (list = ls ())
R.version.string ## "R version 3.6.1 (2019-07-05)"

# libraries
library(dplyr)

# Set working directory
setwd("~/Documents/STRs/ANALYSIS/population_research/MAIN_ANCESTRY/")

# Load unrel genomes batch1
l_unrelated_b1 = read.table("./60k_HWE_30k_random_unrelated_participants.txt", stringsAsFactors = F)
l_unrelated_b1 = l_unrelated_b1$V1
length(l_unrelated_b1)
# 38344

# Load unrel genomes batch2
l_unrelated_b2 = read.table("./batch2/l_unrelated_55603_genomes_batch2.txt", stringsAsFactors = F)
l_unrelated_b2 = l_unrelated_b2$V1
length(l_unrelated_b2)
# 55603

l_unrelated_merged = unique(c(l_unrelated_b1,
                            l_unrelated_b2))
length(l_unrelated_merged)
# 63702

write.table(l_unrelated_merged,
            "./list_63702_UNRELATED_unique_genomes_batch1_batch2.txt",
            quote = F, row.names = F, col.names = F)
