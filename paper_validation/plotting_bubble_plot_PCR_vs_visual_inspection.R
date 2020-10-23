# Objective: plot a bubble plot with the correlation between EH estimations after visual inspection and PCR repeat-sizes
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
setwd("/Users/kibanez/Documents/STRs/VALIDATION/bubble_plots/")

# Load golden validation table - EHv2.5.5
val_data = read.csv("./GEL_accuracy_final_not_NCL.tsv",
                    sep = "\t",
                    header = T,
                    stringsAsFactors = F)
dim(val_data)
# 418  11 - no NCL

output_folder = "./figures/"

group.colors = c("AR" = "#1B9E77", "ATN1" = "#D95F02", "ATXN1" ="#7570B3", "ATXN2" = "#E7298A", "ATXN3" = "#66A61E", 
                 "ATXN7" = "#E6AB02", "CACNA1A" = "#A6761D", "FXN" = "#666666", "HTT" ="#E41A1C", "TBP" = "#FF7F00", 
                 "C9orf72" = "#FFFF33", "FMR1" = "#F781BF", "PPP2R2B" = "black")

# Mike's suggestion - 1
# Truth short and truth long alleles in X axis
# minEhv3 and maxEHv3 in Y axis

# Let's take the important meat
exp_alleles_v2 = c(as.integer(val_data$Truth.Short.Allele), as.integer(val_data$Truth.Long.Allele))
eh_alleles_v2 = c(as.integer(val_data$min.EHv312.a1), as.integer(val_data$max.EHv312.a2))
locus_v2 = c(val_data$locus, val_data$locus)

# Remove NAs
index_NA = which(is.na(exp_alleles_v2))
exp_alleles_v2 = exp_alleles_v2[-index_NA]
eh_alleles_v2 = eh_alleles_v2[-index_NA]
locus_v2 = locus_v2[-index_NA]

# Create dataframe with exp, eh, freq for each locus
df_data_with_freq_v2 = data.frame()
l_locus = unique(locus_v2)
for(i in 1:length(l_locus)){
  aux_validation_a1 = val_data %>% filter(locus %in% l_locus[i]) %>% select(Truth.Short.Allele) %>% pull() %>% as.integer() 
  aux_validation_a2 = val_data %>% filter(locus %in% l_locus[i]) %>% select(Truth.Long.Allele) %>% pull() %>% as.integer() 
  aux_exp_alleles_v2 = c(aux_validation_a1, aux_validation_a2)
  
  aux_eh_a1 = val_data %>% filter(locus %in% l_locus[i]) %>% select(min.EHv312.a1) %>% pull() %>% as.integer() 
  aux_eh_a2 = val_data %>% filter(locus %in% l_locus[i]) %>% select(max.EHv312.a2) %>% pull() %>% as.integer() 
  aux_eh_alleles_v2 = c(aux_eh_a1, aux_eh_a2)
  
  data_aux = xyTable(aux_exp_alleles_v2, aux_eh_alleles_v2)
  
  df_data_aux = data.frame(eh_alleles = data_aux$y,
                           exp_alleles = data_aux$x,
                           number_of_alleles = data_aux$number,
                           locus = rep(l_locus[i], length(data_aux$x)))
  # Concat info per locus
  df_data_with_freq_v2 = rbind(df_data_with_freq_v2,
                               df_data_aux)
  
}

max_value = max(df_data_with_freq_v2$eh_alleles, 
                df_data_with_freq_v2$exp_alleles,
                na.rm = TRUE) + 5
df_strategy1 = df_data_with_freq_v2

# Mike's suggestion - 2
# min and max PCR sizes in X axis
# minEhv3 and maxEHv3 in Y axis

# Let's take the important meat
exp_alleles_v2 = c(as.integer(val_data$exp_PCR_a1), as.integer(val_data$exp_PCR_a2))
eh_alleles_v2 = c(as.integer(val_data$min.EHv312.a1), as.integer(val_data$max.EHv312.a2))
locus_v2 = c(val_data$locus, val_data$locus)

# Remove NAs
index_NA = which(is.na(exp_alleles_v2))
exp_alleles_v2 = exp_alleles_v2[-index_NA]
eh_alleles_v2 = eh_alleles_v2[-index_NA]
locus_v2 = locus_v2[-index_NA]

# Create dataframe with exp, eh, freq for each locus
df_data_with_freq_v2 = data.frame()
l_locus = unique(locus_v2)
for(i in 1:length(l_locus)){
  aux_validation_a1 = val_data %>% filter(locus %in% l_locus[i]) %>% select(exp_PCR_a1) %>% pull() %>% as.integer() 
  aux_validation_a2 = val_data %>% filter(locus %in% l_locus[i]) %>% select(exp_PCR_a2) %>% pull() %>% as.integer() 
  aux_exp_alleles_v2 = c(aux_validation_a1, aux_validation_a2)
  
  aux_eh_a1 = val_data %>% filter(locus %in% l_locus[i]) %>% select(min.EHv312.a1) %>% pull() %>% as.integer() 
  aux_eh_a2 = val_data %>% filter(locus %in% l_locus[i]) %>% select(max.EHv312.a2) %>% pull() %>% as.integer() 
  aux_eh_alleles_v2 = c(aux_eh_a1, aux_eh_a2)
  
  data_aux = xyTable(aux_exp_alleles_v2, aux_eh_alleles_v2)
  
  df_data_aux = data.frame(eh_alleles = data_aux$y,
                           exp_alleles = data_aux$x,
                           number_of_alleles = data_aux$number,
                           locus = rep(l_locus[i], length(data_aux$x)))
  # Concat info per locus
  df_data_with_freq_v2 = rbind(df_data_with_freq_v2,
                               df_data_aux)
  
}

max_value = max(df_data_with_freq_v2$eh_alleles, 
                df_data_with_freq_v2$exp_alleles,
                na.rm = TRUE) + 5

df_strategy2 = df_data_with_freq_v2

group.colors.gray = c("AR" = "#DCDCDC", "ATN1" = "#D3D3D3", "ATXN1" ="#C0C0C0", "ATXN2" = "#A9A9A9", "ATXN3" = "#808080", 
                 "ATXN7" = "#696969", "CACNA1A" = "#778899", "FXN" = "#708090", "HTT" ="#2F4F4F", "TBP" = "#AAAAAA", 
                 "C9orf72" = "#BBBBBB", "FMR1" = "#CCCCCC", "PPP2R2B" = "black")

# Combining joint_plot_mike1 and joint_plot_mike2 into a single one
tontz = ggplot() +
geom_point(data = df_strategy2, aes(x = exp_alleles, y = eh_alleles, size = number_of_alleles), color = "#B8B8B8") +
  geom_point(data = df_strategy1, aes(color = factor(locus), x = exp_alleles, y = eh_alleles, size = number_of_alleles), alpha = 0.7) +  
  xlim(5,max_value) +
  ylim(5,max_value) +
  geom_abline(method = "lm", formula = x ~ y, linetype = 2, colour = "gray") +  
  coord_equal() +
  labs(title = "", 
       y = "EH repeat sizes", 
       x = "PCR repeat sizes") + 
  scale_fill_manual(values=group.colors) +  
  theme_light() +
  theme(legend.title = element_blank(),
        text = element_text(size=13),
        axis.text.x.top = element_text()) +
  guides(size = FALSE)

png("./figures/Figure2B_before_vs_after_visualQC_LANCET_filter_ATN1_ATXN3_from_NCL_600dpi.png",units="in", width=5, height=5, res=600)
print(tontz)
dev.off()

ggsave(file="./figures/Figure2B_before_vs_after_visualQC_LANCET_600dpi.svg", plot=tontz, dpi = 600)


# breakdown by locus
geom_vline(xintercept = l_premut_cutoff[i], colour = 'red', lty = 2) + 
  geom_hline(yintercept = l_premut_cutoff[i], colour = 'red', lty = 2) + 
  
breakdown_by_locus = ggplot(df_strategy1) +
  geom_point(data = df_strategy2, aes(x = exp_alleles, y = eh_alleles, size = number_of_alleles), color = "#B8B8B8") +
  geom_point(data = df_strategy1, aes(color = factor(locus), x = exp_alleles, y = eh_alleles, size = number_of_alleles), alpha = 0.7) +  
  xlim(5,max_value) +
  ylim(5,max_value) +
  geom_abline(method = "lm", formula = x ~ y, linetype = 2, colour = "gray") +  
  coord_equal() +
  labs(title = "", 
       y = "EH repeat sizes", 
       x = "PCR repeat sizes") + 
  scale_fill_manual(values=group.colors) +  
  theme(legend.title = element_blank(),
        axis.text.x.top = element_text()) +
  guides(size = FALSE) +
  facet_wrap(locus~ .)

png("./figures/FigureS3_before_vs_after_visualQC_breakdown_LANCET_filter_ATN1_ATXN3_from_NCL_600dpi.png",units="in", width=5, height=5, res=600)
print(breakdown_by_locus)
dev.off()

ggsave(file="./figures/Figure2_before_vs_after_visualQC_breakdown_LANCET_600dpi.svg", plot=breakdown_by_locus, dpi = 600)
