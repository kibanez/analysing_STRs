# Objetive: represent number of expanded genomes beyond the premutation cut-off together with the PC values
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
setwd("/Users/kibanez/Documents/STRs/ANALYSIS/population_research/PAPER/expanded_genomes_main_pilot/feb2021/beyond_premut/")

# Load data
pc_data = read.csv("~/Documents/STRs/ANALYSIS/population_research/MAIN_ANCESTRY/batch2/aggV2_bedmerge_30KSNPs_labkeyV9_08062020_update_PCsancestryrelated.tsv",
                     stringsAsFactors = F,
                     header = T,
                     sep = " ")
dim(pc_data)
# 78388  17

popu_data = read.csv("~/Documents/STRs/ANALYSIS/population_research/MAIN_ANCESTRY/batch2/aggV2_M30K_60K_1KGP3_ancestry_assignment_probs_R9_08062020.tsv",
                     stringsAsFactors = F,
                     header = T,
                     sep = " ")
dim(popu_data)
# 78388 33

popu_merged = left_join(pc_data,
                        popu_data %>% select(plate_key, ancestry0_8),
                        by = "plate_key")
dim(popu_merged)
# 78388  18

# Load the final table with the number of expanded genomes beyond the premutation threshold
table1 = read.csv("./13loci_beyond_premut_cutoff_to_review_VGD_enriched_pathoFinalDecision_100621.tsv",
                  stringsAsFactors = F,
                  sep = "\t",
                  header = T)
dim(table1)
# 2696  9

# Merge with PC values
table1 = left_join(table1,
                   popu_merged %>% select(plate_key, Pc1, Pc2, Pc3, Pc4, Pc5, Pc6, Pc7, Pc8, Pc9, Pc10),
                   by = c("platekey" = "plate_key"))

l_genes = unique(table1$locus)

table1_locus = table1 %>%
  filter(locus %in% l_genes[i], is_unrel, is_125 %in% "No", Final.decision %in% "Yes")

# Let's plot
ggplot(data=popu_batch2, 
       aes(x=PC2, y=PC1, colour = ancestry0_8)) +
  geom_hex(bins=300) +
  xlab("PC2 across 78,388 genomes") +
  ylab("PC1 across 78,388 genomes") +
  guides(fill = FALSE)



