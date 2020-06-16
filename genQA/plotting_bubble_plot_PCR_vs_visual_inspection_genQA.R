# Objective: plot a bubble plot with the correlation between EH and genQA original labo repeat-sizes
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
setwd("/Users/kibanez/Documents/STRs/VALIDATION/genQA/genQA/")

# Load golden validation table - EHv2.5.5
val_data = read.csv("./merged_batch1_batch3_STRs.tsv",
                    sep = "\t",
                    header = T,
                    stringsAsFactors = F)
dim(val_data)
# 46  17

# 1 - Define min and max alleles for both PCR and EH 
# Don't work for AR, FMR1 - when we do have NA values
val_data = val_data %>% 
  group_by(Platekey, locus) %>%
  mutate(min_PCR = min(GenQA_a1, GenQA_a2, na.rm = F),
         max_PCR = max(GenQA_a1, GenQA_a2, na.rm = F),
         min_EH = min(GEL_a1, GEL_a2, na.rm = F),
         max_EH = max(GEL_a1, GEL_a2, na.rm = F)) %>%
  ungroup() %>%
  as.data.frame()
dim(val_data)
# 46  21

# Let's take the important meat: experimentally validated data and EH estimations
exp_alleles_v2 = c(as.integer(val_data$min_PCR), as.integer(val_data$max_PCR))
eh_alleles_v2 = c(as.integer(val_data$min_EH), as.integer(val_data$max_EH))
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
  aux_validation_a1 = val_data %>% filter(locus %in% l_locus[i]) %>% select(GenQA_a1) %>% pull() %>% as.integer() 
  aux_validation_a2 = val_data %>% filter(locus %in% l_locus[i]) %>% select(GenQA_a2) %>% pull() %>% as.integer() 
  aux_exp_alleles_v2 = c(aux_validation_a1, aux_validation_a2)
  
  aux_eh_a1 = val_data %>% filter(locus %in% l_locus[i]) %>% select(GEL_a1) %>% pull() %>% as.integer() 
  aux_eh_a2 = val_data %>% filter(locus %in% l_locus[i]) %>% select(GEL_a2) %>% pull() %>% as.integer() 
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

output_folder = "./figures/"

max_value = max(df_data_with_freq_v2$eh_alleles, 
                df_data_with_freq_v2$exp_alleles,
                na.rm = TRUE) + 5

group.colors = c("AR" = "#1B9E77", "ATN1" = "#D95F02", "ATXN1" ="#7570B3", "ATXN2" = "#E7298A", "ATXN3" = "#66A61E", 
                 "ATXN7" = "#E6AB02", "CACNA1A" = "#A6761D", "FXN" = "#666666", "HTT" ="#E41A1C", "TBP" = "#FF7F00", 
                 "C9orf72" = "#FFFF33", "FMR1" = "#F781BF", "PPP2R2B" = "black")

# Sort by locus
df_data_with_freq_v2 = df_data_with_freq_v2[order(df_data_with_freq_v2$locus),]

joint_plot = ggplot() +
geom_point(data = df_data_with_freq_v2, aes(x = exp_alleles, y = eh_alleles, size = number_of_alleles, color = factor(locus))) +
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

png("./figures/genQA_PCR_vs_EH_generalView.png",units="in", width=5, height=5, res=300)
print(joint_plot)
dev.off()

# breakdown by locus
breakdown_by_locus = ggplot() +
  geom_point(data = df_data_with_freq_v2, aes(x = exp_alleles, y = eh_alleles, size = number_of_alleles, color = factor(locus))) +
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
  guides(size = FALSE) +
  facet_wrap(locus~ .)

png("./figures/genQA_PCR_vs_EH_generalView_breakdown.png",units="in", width=5, height=5, res=600)
print(breakdown_by_locus)
dev.off()

# Breakdown by locus and cutoff
# For each locus
l_genes = c("AR", "ATN1", "ATXN1", "ATXN2", "ATXN7", "CACNA1A", "C9orf72", "DMPK", "HTT", "FMR1", "FXN", "TBP", "PPP2R2B")
l_premut_cutoff = c(34,34,35,31,17,17,30,50,35,55,44,41,51)


for(i in 1:length(l_genes)){
  
  df_data_with_freq_v2_indiv = df_data_with_freq_v2 %>% 
    filter(locus %in% l_genes[i])
  
  max_value_indiv = max(df_data_with_freq_v2_indiv$eh_alleles, 
                        df_data_with_freq_v2_indiv$exp_alleles,
                        na.rm = TRUE) + 5
  
  joint_plot_individual = ggplot() +
    geom_point(data = df_data_with_freq_v2_indiv,
               aes(x = exp_alleles, y = eh_alleles, size = number_of_alleles)) +
    #geom_rect(aes(xmin=l_premut_cutoff[i], xmax=max_value_indiv, ymin=5,ymax=max_value_indiv), alpha=0.2, fill="red") +
    labs(title = l_genes[i], 
         y = "EH repeat sizes", 
         x = "PCR repeat sizes") + 
    scale_fill_manual(values=group.colors.classi) +  
    theme_light() +
    theme(legend.title = element_blank(),
          text = element_text(size=13),
          axis.text.x.top = element_text(),
          aspect.ratio=1) +
    geom_vline(xintercept = l_premut_cutoff[i], colour = 'red', lty = 2) + 
    geom_hline(yintercept = l_premut_cutoff[i], colour = 'red', lty = 2) + 
    guides(size = FALSE) + 
    coord_equal() +
    xlim(5,max_value_indiv) +
    ylim(5,max_value_indiv) +
    geom_abline(method = "lm", formula = x ~ y, linetype = 2, colour = "gray") 
  
  
  png(paste(paste("./figures/genQA_PCR_vs_EH_individual", l_genes[i], sep = "_"), "" , sep = ".png"),units="in", width=5, height=5, res=300)
  print(joint_plot_individual)
  dev.off()
  
}


# Now with classification
# Enrich with classification
df_classi = data.frame()
for (i in 1:length(l_genes)){
  aux = df_data_with_freq_v2 %>% filter(locus %in% l_genes[i])
  
  aux_classi = ifelse(aux$exp_alleles > l_premut_cutoff[i], "TP", "TN")
  
  aux$classi = aux_classi
  df_classi = rbind(df_classi,
                    aux)
}

group.colors.classi = c("TN" = "#12C420","TP" = "#C43212")

for(i in 1:length(l_genes)){
  
  df_classi_indiv = df_classi %>% 
    filter(locus %in% l_genes[i], !is.na(classi))
  
  max_value_indiv = max(df_classi_indiv$eh_alleles, 
                        df_classi_indiv$exp_alleles,
                        na.rm = TRUE) + 5
  
  joint_plot_individual = ggplot() +
    geom_point(data = df_classi_indiv,
               aes(x = exp_alleles, y = eh_alleles, size = number_of_alleles, color = factor(classi))) +
    #geom_rect(aes(xmin=l_premut_cutoff[i], xmax=max_value_indiv, ymin=5,ymax=max_value_indiv), alpha=0.2, fill="red") +
    labs(title = l_genes[i], 
         y = "EH repeat sizes", 
         x = "PCR repeat sizes") + 
    scale_fill_manual(values=group.colors.classi) +  
    theme_light() +
    theme(legend.title = element_blank(),
          text = element_text(size=13),
          axis.text.x.top = element_text()) +
    geom_vline(xintercept = l_premut_cutoff[i], colour = 'red', lty = 2) + 
    geom_hline(yintercept = l_premut_cutoff[i], colour = 'red', lty = 2) + 
    guides(size = FALSE) + 
    coord_equal() +
    xlim(5,max_value_indiv) +
    ylim(5,max_value_indiv) +
    geom_abline(method = "lm", formula = x ~ y, linetype = 2, colour = "gray") 
  
  
  png(paste(paste("./figures/genQA_PCR_vs_EH_individual_with_classi", l_genes[i], sep = "_"), "" , sep = ".png"),units="in", width=5, height=5, res=300)
  print(joint_plot_individual)
  dev.off()
  
}

