# Objective: analyse AR epi across EUR genomes in the 1Kg project
date ()
Sys.info ()[c("nodename", "user")]
commandArgs ()
rm (list = ls ())
R.version.string ## "R version 3.6.1 (2019-07-05)"

library(ggplot2); packageDescription ("ggplot2", fields = "Version") #"3.2.1"
library(grid); packageDescription ("grid", fields = "Version") #"3.6.2"
library(gridExtra); packageDescription ("gridExtra", fields = "Version") #"2.3"
library(reshape2); packageDescription ("reshape2", fields = "Version") #"1.4.3"
require(dplyr); packageDescription ("dplyr", fields = "Version") #"0.8.3"

# Set working dir
setwd("~/Documents/STRs/ANALYSIS/AR_kennedy/1Kg/")

# Load data
all_merged = read.csv("~/Documents/STRs/ANALYSIS/population_research/1kg/data/all_data/merged/merged_all_2504_genomes_1Kg_EHv3.2.2.tsv",
                      stringsAsFactors = F,
                      header = T,
                      sep = "\t")
dim(all_merged)
# 1342  12

# Let's select `AR` locus
ar_merged = all_merged %>%
  filter(gene %in% "AR") %>%
  select(gene, allele, num_samples)
ar_merged = unique(ar_merged)
dim(ar_merged)
# 30  3

# Interesting thing: there are 4 genomes with an expanded AR expansion
ar_merged %>% filter(allele > 35)
#gene allele num_samples
#1   AR     36           1
#2   AR     39           1
#3   AR     38           1
#4   AR     37           1
 


# Plot percentage of each allele within each subpopulation
violin_plot = ggplot(eur_merged, aes(x = subpopu, y=allele, fill = subpopu)) +
  geom_violin() +
  coord_flip() +
  xlab("Sub-population across EUR in 1Kg dataset") + 
  ylab("Repeat sizes (repeat units)") + 
  ggtitle("AR across EUR sub-populations within 1Kg") +
  geom_boxplot(width=0.1) +
  guides(fill=guide_legend(reverse=TRUE)) 

violin_output_png = "./plots/violin_plots_across_EUR_subpopus.png"
png(violin_output_png, units="in", width=5, height=5, res=300)
print(violin_plot)
dev.off()

# Plot joint ancestry repeat-size distributions
min_value = min(eur_merged$allele, 34)
max_value = max(eur_merged$allele, 38)
joint_plot = ggplot(eur_merged, aes(x = allele, y = percent_subpopu, group = subpopu, color = subpopu)) + 
  geom_line() + 
  ylab("Allele frequency") + 
  xlab("Repeat sizes (repeat units)") + 
  geom_vline(xintercept = 34, colour = 'blue', lty = 2) + 
  geom_vline(xintercept = 38, colour = 'red', lty = 2) + 
  coord_cartesian(xlim = c(min_value,max_value))

joint_ancestries = "./plots/joint_subpopus_in_AR_1Kg.png"
png(joint_ancestries, units="in", width=5, height=5, res=300)
print(joint_plot)
dev.off()
