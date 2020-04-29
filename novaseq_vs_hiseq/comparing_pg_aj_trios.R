# Objective: plot a bubble plot or 1:1 plot to analyse Hiseq vs Novaseq performance from STRs point of view
date ()
Sys.info ()[c("nodename", "user")]
commandArgs ()
rm (list = ls ())
R.version.string ## "R version 3.6.3(2020-02-29)"

library(ggplot2); packageDescription ("ggplot2", fields = "Version") #"2.2.1"
library(grid); packageDescription ("grid", fields = "Version") #"3.3.2"
library(gridExtra); packageDescription ("gridExtra", fields = "Version") #"2.3"
library(reshape2); packageDescription ("reshape2", fields = "Version") #"1.4.3"
require(dplyr); packageDescription ("dplyr", fields = "Version") #"0.7.4"
library(lattice); packageDescription ("lattice", fields = "Version") #"0.20-35"
library(tidyverse); packageDescription ("tidyverse", fields = "Version") #"1.2.1"
library(RColorBrewer); packageDescription ("RColorBrewer", fields = "Version") #"1.1-2"
library(cowplot); packageDescription ("cowplot", fields = "Version") #"1.0.0"

# Set working environment
setwd("/Users/kibanez/Documents/STRs/ANALYSIS/novaseq/STRs/")

# Load data

# Hiseq data
ash_hiseq = read.csv("./HG002_Edico.HiSeq.repeats.tsv",
                     stringsAsFactors = F, 
                     sep = "\t",
                     header = F)
ash_hiseq$tech = rep("HiSeq", length(ash_hiseq$V1))
ash_hiseq$trio = rep("Ash", length(ash_hiseq$V1))

pg_hiseq = read.csv("./NA12878_S1_Edico.HiSeq.repeats.tsv",
                    stringsAsFactors = F,
                    sep = "\t",
                    header = F)
pg_hiseq$tech = rep("HiSeq", length(pg_hiseq$V1))
pg_hiseq$trio = rep("PG", length(pg_hiseq$V1))


# Novaseq1
ash_nova1 = read.csv("./LP4100018-DNA_E11_Edico.Ash.NovaSeq.repeats.tsv",
                     stringsAsFactors = F,
                     sep = "\t",
                     header = F)
ash_nova1$tech = rep("NovaSeq", length(ash_nova1$V1))
ash_nova1$trio = rep("Ash", length(ash_nova1$V1))

pg_nova1= read.csv("./LP4100018-DNA_E05_Edico.NovaSeq.repeats.tsv",
                    stringsAsFactors = F, 
                    sep = "\t",
                    header = F)
pg_nova1$tech = rep("NovaSeq", length(pg_nova1$V1))
pg_nova1$trio = rep("PG", length(pg_nova1$V1))

# Novaseq2
ash_nova2 = read.csv("./LP4100018-DNA_F11_Edico.Ash.NovaSeq.repeats.tsv",
                     stringsAsFactors = F,
                     sep = "\t",
                     header = F)
ash_nova2$tech = rep("NovaSeq", length(ash_nova2$V1))
ash_nova2$trio = rep("Ash", length(ash_nova2$V1))

pg_nova2= read.csv("./LP4100016-DNA_D02_Edico.NovaSeq.repeats.tsv",
                    stringsAsFactors = F, 
                    sep = "\t",
                    header = F)
pg_nova2$tech = rep("NovaSeq", length(pg_nova2$V1))
pg_nova2$trio = rep("Ash", length(pg_nova2$V1))


merged_all = rbind(ash_hiseq,
                   pg_hiseq,
                   ash_nova1,
                   pg_nova1,
                   ash_nova2,
                   pg_nova2)
colnames(merged_all) = c("sequencing", "locus", "EHv3_a1", "EHv3_a2", "technology", "trio")

# Focus only in the 13 loci we have validated clinically (paper)
l_genes = c("AR_CAG", "ATN1_CAG", "ATXN1_CAG", "ATXN2_CAG", "ATXN3_CAG", "ATXN7_CAG", "CACNA1A_CAG", "FXN_GAA", "HTT_CAG", "C9orf72_GGGGCC", "FXN_GAA","TBP_CAG", "PPPP2R2B_CAG")

# Filling manually colours (locus)
brewer.pal(n = 8, name = "Dark2")
brewer.pal(n = 8, name = "Set1")

group.colors = c("AR_CAG" = "#1B9E77", "ATN1_CAG" = "#D95F02", "ATXN1_CAG" ="#7570B3", "ATXN2_CAG" = "#E7298A", "ATXN3_CAG" = "#66A61E", 
                 "ATXN7_CAG" = "#E6AB02", "CACNA1A_CAG" = "#A6761D", "FXN_GAA" = "#666666", "HTT_CAG" ="#E41A1C", "TBP_CAG" = "#FF7F00", 
                 "C9orf72_GGGGCC" = "#FFFF33", "PPP2R2B_CAG" = "black")

merged_genes = merged_all %>%
  filter(locus %in% l_genes)

ash_hiseq = merged_genes %>%
  filter(trio %in% "Ash", technology %in% "HiSeq")
ash_novaseq = merged_genes %>%
  filter(trio %in% "Ash", technology %in% "NovaSeq")

pg_hiseq = merged_genes %>%
  filter(trio %in% "PG", technology %in% "HiSeq")

joint_plot = ggplot(merged_genes, 
                    aes(x = ash_hiseq, y = ash_novaseq), colour = factor(locus)) + 
  geom_point(aes(fill = factor(locus))) + 
  labs(title = "", 
       y = "Repeat sizes for each allele \n NovaSeq", 
       x = "Repeat sizes for each allele \n HiSeq") + 
  geom_abline(method = "lm", formula = x ~ y, linetype = 2, colour = "gray") +  
  coord_equal() +
  scale_fill_manual(values=group.colors) +  
  theme(legend.title = element_blank(),
        axis.text.x.top = element_text()) + 
  guides(size = FALSE) 

