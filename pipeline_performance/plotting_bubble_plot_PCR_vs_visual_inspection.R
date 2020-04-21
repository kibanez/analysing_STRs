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
val_data = read.csv("EHv2_avg_VS_EHv2_maxCI/STRVALIDATION_ALLDATA_2019-10-7_ALL_kibanez_EHv255_avg_VS_EHv255_maxCI_checkFXN_withPileup_and_expValidatedData_tweaking_ATN1_updated_AR_from_NHNN.tsv",
                    sep = "\t",
                    header = T,
                    stringsAsFactors = F)
dim(val_data)
# 635  22

# Filter the good ones
# 1 - Only keep with `Pileup_quality` == good or Good, `MISSING` and `blanks`
val_data = val_data %>%
  filter(Pileup_quality %in% "Good" | Pileup_quality %in% "good" | Pileup_quality %in% "" | Pileup_quality %in% "MISSING")
dim(val_data)
# 543  22

# Let's simplify the data we need from `val_data`
val_data = val_data %>%
  select(LP_Number, locus_bioinfo, locus, STR_a1, STR_a2, EH_a1_avg, EH_a2_avg)

# 2 - Only keep experimental val numbers (STR_a1, STR_a2) that are integer
# Let's see index for which `STR_a1` and `STR_a2` separately have not integer or numbers for allele estimation
index_a1 = c(which(val_data$STR_a1 == "normal"), 
             which(val_data$STR_a1 == "positive"), 
             which(is.na(val_data$STR_a1)))
length(index_a1)
# 8

index_a2 = c(which(val_data$STR_a2 == "positive"), 
             which(val_data$STR_a2 == "premutation"), 
             which(val_data$STR_a2 == "na"), 
             which(val_data$STR_a2 == "full_mutation"), 
             which(val_data$STR_a2 == "EXP"),
             which(val_data$STR_a2 == "."),
             which(is.na(val_data$STR_a2)))
length(index_a2)
# 43

# EH alleles with no integer
index_eh_a1 = c(which(is.na(val_data$EH_a1_avg)))
length(index_eh_a1)
# 1

index_eh_a2 = c(which(is.na(val_data$EH_a2_avg)))
length(index_eh_a2)
# 12

# We will use `index_a1` and `index_a2` when filtering out STR_a1, and EH_a1_avg; and STR_a2 and EH_a2_avg respectively
index_allele1 = unique(c(index_eh_a1, index_a1))
length(index_allele1)
# 9

index_allele2 = unique(c(index_eh_a2, index_a2))
length(index_allele2)
# 45

# in case there is only 1 allele, we should have 1 allele for expValidation and EH
# in case there are 2 alleles, then the `a1` for expValidation and EH should be the minimum (minor repeat size)
val_data$STR_a1 = as.integer(val_data$STR_a1)
val_data$STR_a2 = as.integer(val_data$STR_a2)
val_data$EH_a1_avg = as.integer(val_data$EH_a1_avg)
val_data$EH_a2_avg = as.integer(val_data$EH_a2_avg)
for (i in 1:length(val_data$LP_Number)){
  val_validation_a1 = val_data$STR_a1[i]
  val_validation_a2 = val_data$STR_a2[i]
  val_eh_a1 = val_data$EH_a1_avg[i]
  val_eh_a2 = val_data$EH_a2_avg[i]
  
  # For alleles that there is no estimation (EH only calls 1 allele) we do have a `0` -- we ignore/avoid this step in these cases
  if (val_eh_a2 == 0 | is.na(val_eh_a2)){
    min_validation = val_validation_a1
    min_eh = val_eh_a1
    max_validation = NA
    max_eh = NA
  }else{
    min_validation = min(val_validation_a1, val_validation_a2)
    max_validation = max(val_validation_a1, val_validation_a2)
    min_eh = min(val_eh_a1, val_eh_a2)
    max_eh = max(val_eh_a1, val_eh_a2)
  }
  # Post-processing the new values
  val_data$validation_a1[i] = min_validation
  val_data$validation_a2[i] = max_validation
  val_data$eh_a1[i] = min_eh
  val_data$eh_a2[i] = max_eh
}
dim(val_data)
# 543  11

# Let's take the important meat: experimentally validated data and EH estimations
exp_alleles_v2 = c(as.integer(val_data$validation_a1), as.integer(val_data$validation_a2))
eh_alleles_v2 = c(as.integer(val_data$eh_a1), as.integer(val_data$eh_a2))
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
  aux_validation_a1 = val_data %>% filter(locus %in% l_locus[i]) %>% select(validation_a1) %>% pull() %>% as.integer() 
  aux_validation_a2 = val_data %>% filter(locus %in% l_locus[i]) %>% select(validation_a2) %>% pull() %>% as.integer() 
  aux_exp_alleles_v2 = c(aux_validation_a1, aux_validation_a2)
  
  aux_eh_a1 = val_data %>% filter(locus %in% l_locus[i]) %>% select(eh_a1) %>% pull() %>% as.integer() 
  aux_eh_a2 = val_data %>% filter(locus %in% l_locus[i]) %>% select(eh_a2) %>% pull() %>% as.integer() 
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
       y = "Repeat sizes for each allele \n Expansion Hunter (EH-v2.5.5)", 
       x = "Repeat sizes for each allele \n Experimental validation") + 
  geom_abline(method = "lm", formula = x ~ y, linetype = 2, colour = "gray") +  
  coord_equal() +
  scale_fill_manual(values=group.colors) +  
  theme(legend.title = element_blank(),
        axis.text.x.top = element_text()) + 
  guides(size = FALSE) 


png("figures/joint_bubble_plot_EHv2_generalView.png", units="in", width=5, height=5, res=300)
print(joint_plot)
dev.off()

pdf("figures/joint_bubble_plot_EHv2_generalView.pdf")
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

