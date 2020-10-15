# Libraries
library(dplyr)
library(ggplot2)
#install.packages("svglite")

# Set working environment
setwd("~/Documents/STRs/PAPERS/VALIDATION_PAPER/LANCET/figures/Figure3/")

# Load data
mega_merge_no_NA = read.csv("./mega_merge_no_NA.tsv",
                            stringsAsFactors = F,
                            header = T,
                            sep = "\t")
dim(mega_merge_no_NA)
# 106004  19

# Reformat gene names
mega_merge_no_NA$gene = gsub("_CAG", "", mega_merge_no_NA$gene)

# add normal and pathogenic thresholds
patho_thresholds = data.frame(gene = c("ATN1", "ATXN2", "ATXN7", "HTT"),
                              patho_cutoff2 = c(48, 33, 36, 40))

mega_merge_no_NA = left_join(mega_merge_no_NA,
                             patho_thresholds,
                             by = "gene")


vline.data = as_tibble(patho_thresholds)

fig3 = ggplot(mega_merge_no_NA, aes( x= repeat_size), fill= gene) + 
  geom_histogram(position="identity", binwidth = 1, size = 0.4, colour = "black", fill="skyblue2") + 
  scale_x_continuous(name = "Repeat size", breaks = seq(0,130,10), expand = c(0, 0)) + 
  scale_y_continuous(name = "Allele Count", expand = c(0, 0)) + 
  coord_cartesian(ylim=c(0, 3500))  + 
  theme_light() + theme(text = element_text(size=16)) + 
  facet_grid(gene ~ .) +
  geom_vline(aes(xintercept = patho_cutoff2), vline.data, colour = "red", linetype="dotted")
  

ggsave(file="Figure3_600dpi.svg", plot=fig3, dpi = 600)

