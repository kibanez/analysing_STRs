# Objetive: analyse population release - batch2
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

# Load popu data with related/no related info
popu_batch2 = read.csv("./MAIN_ANCESTRY/batch2/aggV2_bedmerge_30KSNPs_labkeyV9_08062020_update_PCsancestryrelated.tsv",
                       sep = " ",
                       stringsAsFactors = F,
                       header = T)
dim(popu_batch2)
# 78388  17

unique(length(popu_batch2$plate_key))
# 78388

# List of unrelated genomes
l_unrelated = popu_batch2 %>% filter(unrelated_set == 1) %>% select(plate_key) %>% unique() %>% pull()
length(l_unrelated)
# 55847

write.table(l_unrelated,
            "./MAIN_ANCESTRY/batch2/l_unrelated_55847_genomes_batch2.txt",
            quote = F,
            row.names = F,
            col.names = F)

