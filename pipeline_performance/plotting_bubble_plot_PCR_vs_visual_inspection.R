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

# Set working environment
setwd("/Users/kibanez/Documents/STRs/VALIDATION/bubble_plots/")

# Load golden validation table - EHv2.5.5
val_data = read.csv("./GEL_accuracy_final.tsv",
                    sep = "\t",
                    header = T,
                    stringsAsFactors = F)
dim(val_data)
# 616  9

# 1 - Define min and max alleles for both PCR and EH 
val_data = val_data %>% 
  group_by(LP_number, locus) %>%
  mutate(min_PCR = min(exp_PCR_a1, exp_PCR_a2),
         max_PCR = max(exp_PCR_a1, exp_PCR_a2),
         min_EH = min(Truth.Short.Allele, Truth.Long.Allele),
         max_EH = max(Truth.Short.Allele, Truth.Long.Allele)) %>%
  ungroup() %>%
  as.data.frame()
dim(val_data)
# 616  13

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
  aux_validation_a1 = val_data %>% filter(locus %in% l_locus[i]) %>% select(min_PCR) %>% pull() %>% as.integer() 
  aux_validation_a2 = val_data %>% filter(locus %in% l_locus[i]) %>% select(max_PCR) %>% pull() %>% as.integer() 
  aux_exp_alleles_v2 = c(aux_validation_a1, aux_validation_a2)
  
  aux_eh_a1 = val_data %>% filter(locus %in% l_locus[i]) %>% select(min_EH) %>% pull() %>% as.integer() 
  aux_eh_a2 = val_data %>% filter(locus %in% l_locus[i]) %>% select(max_EH) %>% pull() %>% as.integer() 
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

#df_data_with_freq_v2$exp_alleles[133] = 100

max_value = max(df_data_with_freq_v2$eh_alleles, 
                df_data_with_freq_v2$exp_alleles,
                na.rm = TRUE) + 5

# Filling manually colours (locus)
brewer.pal(n = 8, name = "Dark2")
# "#1B9E77" "#D95F02" "#7570B3" "#E7298A" "#66A61E" "#E6AB02" "#A6761D" "#666666"
brewer.pal(n = 8, name = "Set1")
# "#E41A1C" "#377EB8" "#4DAF4A" "#984EA3" "#FF7F00" "#FFFF33" "#A65628" "#F781BF"

group.colors = c("AR" = "#1B9E77", "ATN1" = "#D95F02", "ATXN1" ="#7570B3", "ATXN2" = "#E7298A", "ATXN3" = "#66A61E", 
                 "ATXN7" = "#E6AB02", "CACNA1A" = "#A6761D", "FXN" = "#666666", "HTT" ="#E41A1C", "TBP" = "#FF7F00", 
                 "C9orf72" = "#FFFF33", "FMR1" = "#F781BF", "PPP2R2B" = "black")

joint_plot = ggplot(df_data_with_freq_v2, 
                    aes(x = exp_alleles, y = eh_alleles, colour = factor(locus))) + 
  geom_point(aes(fill = factor(locus), size = number_of_alleles)) + 
  xlim(5,max_value) + 
  ylim(5,max_value) + 
  labs(title = "", 
       y = "Repeat sizes for each allele \n Expansion Hunter after visual inspection", 
       x = "Repeat sizes for each allele \n Experimental PCR validation") + 
  geom_abline(method = "lm", formula = x ~ y, linetype = 2, colour = "gray") +  
  coord_equal() +
  scale_fill_manual(values=group.colors) +  
  theme(legend.title = element_blank(),
        axis.text.x.top = element_text()) + 
  guides(size = FALSE) 


png("figures/joint_bubble_plot_EH_visualInspection_vs_PCR_generalView.png", units="in", width=5, height=5, res=300)
print(joint_plot)
dev.off()

pdf("figures/joint_bubble_plot_EH_visualInspection_vs_PCR_generalView.pdf")
print(joint_plot)
dev.off()

# Let's focus on those expansions >70 and <70
# Larger than 70
df_data_with_freq_v2_larger70 = df_data_with_freq_v2 %>% filter(eh_alleles > 70)
max_value = max(df_data_with_freq_v2_larger70$eh_alleles, 
                df_data_with_freq_v2_larger70$exp_alleles,
                na.rm = TRUE) + 5
min_value = min(df_data_with_freq_v2_larger70$eh_alleles, 
                df_data_with_freq_v2_larger70$exp_alleles,
                na.rm = TRUE) - 5

joint_plot_larger70 = ggplot(df_data_with_freq_v2_larger70, 
                    aes(x = exp_alleles, y = eh_alleles, colour = factor(locus))) + 
  geom_point(aes(fill = factor(locus))) + 
  xlim(min_value - 5 ,max_value + 5) + 
  ylim(min_value - 5 ,max_value + 5) + 
  labs(title = "", 
       y = "Repeat sizes for each allele - larger 70 repeats \n Expansion Hunter (EH-v2.5.5)", 
       x = "Repeat sizes for each allele \n Experimental validation") + 
  geom_abline(method = "lm", formula = x ~ y, linetype = 2, colour = "gray") +  
  coord_equal() +
  scale_fill_manual(values=group.colors) +  
  theme(legend.title = element_blank()) + 
  guides(size = FALSE)

png("figures/joint_bubble_plot_EHv2_larger70.png", units="in", width=5, height=5, res=300)
print(joint_plot_larger70)
dev.off()
pdf("figures/joint_bubble_plot_EHv2_larger70.pdf")
print(joint_plot_larger70)
dev.off()

# smaller than 70
df_data_with_freq_v2_smaller70 = df_data_with_freq_v2 %>% filter(eh_alleles <= 70)
max_value = max(df_data_with_freq_v2_smaller70$eh_alleles, 
                df_data_with_freq_v2_smaller70$exp_alleles,
                na.rm = TRUE) + 5
min_value = min(df_data_with_freq_v2_smaller70$eh_alleles, 
                df_data_with_freq_v2_smaller70$exp_alleles,
                na.rm = TRUE) - 5

joint_plot_smaller70 = ggplot(df_data_with_freq_v2_smaller70, 
                             aes(x = exp_alleles, y = eh_alleles, colour = factor(locus))) + 
  geom_point(aes(fill = factor(locus), size = number_of_alleles)) + 
  xlim(5,max_value) + 
  ylim(5,max_value) + 
  labs(title = "", 
       y = "Repeat sizes for each allele - smaller or eq to 70 repeats \n Expansion Hunter (EH-v2.5.5)", 
       x = "Repeat sizes for each allele \n Experimental validation") + 
  geom_abline(method = "lm", formula = x ~ y, linetype = 2, colour = "gray") +  
  coord_equal() +
  scale_fill_manual(values=group.colors) +  
  theme(legend.title = element_blank()) + 
  guides(size = FALSE)

png("figures/joint_bubble_plot_EHv2_smaller70.png", units="in", width=5, height=5, res=300)
print(joint_plot_smaller70)
dev.off()
pdf("figures/joint_bubble_plot_EHv2_smaller70.pdf")
print(joint_plot_smaller70)
dev.off()



# let's plot each locus independently
for(i in 1:length(l_locus)){
  df_data_with_freq_v2_locus = df_data_with_freq_v2 %>% 
    filter(locus %in% l_locus[i])
  
  max_value = max(df_data_with_freq_v2_locus$eh_alleles, 
                  df_data_with_freq_v2_locus$exp_alleles) + 5
  
  file_name = paste(l_locus[i], "experimental_vs_EHv2.5.5", sep = "_")
  pdf_name = paste(file_name, "pdf", sep = ".")
  png_name = paste(file_name, "png", sep = ".")
  pdf_output = paste(output_folder, pdf_name, sep = "")
  png_output = paste(output_folder, png_name, sep = "")
  
  locus_bubble = ggplot(df_data_with_freq_v2_locus, 
         aes(x = exp_alleles, y = eh_alleles)) + 
    geom_point(aes(color = locus, size = number_of_alleles), show.legend = FALSE) + 
    xlim(5,max_value) + 
    ylim(5,max_value) + 
    labs(title = l_locus[i], 
         y = "Repeat sizes for each allele \n Expansion Hunter (EH-v2.5.5)", 
         x = "Repeat sizes for each allele \n Experimental validation") + 
    theme(plot.title = element_text(hjust = 0.5)) +
    geom_abline(method = "lm", formula = x ~ y, linetype = 2, colour = "gray") +  
    coord_equal()
  
  pdf(pdf_output)
  print(locus_bubble)
  dev.off()
  
  png(png_output, units="in", width=5, height=5, res=300)
  print(locus_bubble)
  dev.off()
}

