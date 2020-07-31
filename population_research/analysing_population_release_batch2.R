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
setwd("/Users/kibanez/Documents/STRs/ANALYSIS/population_research/MAIN_ANCESTRY/batch2/")

# Load popu data with related/no related info
popu_batch2 = read.csv("aggV2_bedmerge_30KSNPs_labkeyV9_08062020_update_PCsancestryrelated.tsv",
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
            "l_unrelated_55847_genomes_batch2.txt",
            quote = F,
            row.names = F,
            col.names = F)

# Load popu data with ancestry info
popu_batch2 = read.csv("aggV2_M30K_60K_1KGP3_ancestry_assignment_probs_R9_08062020.tsv",
                       sep = " ",
                       stringsAsFactors = F,
                       header = T)
dim(popu_batch2)
# 78388  33

# Let's plot raw ancestry data
raw_numbers_popus = as.data.frame(table(popu_batch2$ancestry0_8))
colnames(raw_numbers_popus) = c("population", "Number of genomes")


png("figures/barplot_pure_ancestry_groups_raw_numbers.png")
ggplot(raw_numbers_popus, 
       aes(x = reorder(population, -`Number of genomes`), y = `Number of genomes`)) + 
  geom_bar(stat = "identity", aes(fill = population)) + 
  geom_text(aes(label=`Number of genomes`), vjust=-0.5, size = 4, colour = "grey") +
  ylab("Number of genomes") + 
  xlab("Ancestry cohorts - EHv3.2.2 - 78,388 total genomes") 
dev.off()

# Load popu data with PCs
popu_batch2_pc = read.csv("aggV2_M30K_60K_1KGP3projectionPCs_R9_08062020.tsv",
                          sep = " ",
                          stringsAsFactors = F,
                          header = T)
dim(popu_batch2_pc)
# 78388  21


popu_batch2 = left_join(popu_batch2,
                        popu_batch2_pc,
                        by = "plate_key")
dim(popu_batch2)
# 78388  53

png("figures/population_distribution_78388_all_ancestries.png")
ggplot(data=popu_batch2, 
       aes(x=PC2, y=PC1, colour = ancestry0_8)) +
  geom_hex(bins=300) +
  xlab("PC2 across 78,388 genomes") +
  ylab("PC1 across 78,388 genomes") +
  guides(fill = FALSE)
dev.off()
