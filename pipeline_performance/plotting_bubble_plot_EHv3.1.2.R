# Objective: plot a bubble plot with the correlation between EHv2/EHv3 estimations and the experimental validation
date ()
Sys.info ()[c("nodename", "user")]
commandArgs ()
rm (list = ls ())
R.version.string ## "R version 3.6.1 (2019-07-05)"

library(ggplot2); packageDescription ("ggplot2", fields = "Version") #"2.2.1"
library(grid); packageDescription ("grid", fields = "Version") #"3.3.2"
library(gridExtra); packageDescription ("gridExtra", fields = "Version") #"2.3"
library(reshape2); packageDescription ("reshape2", fields = "Version") #"1.4.3"
require(dplyr); packageDescription ("dplyr", fields = "Version") #"0.7.4"
library(lattice); packageDescription ("lattice", fields = "Version") #"0.20-35"
library(tidyverse); packageDescription ("tidyverse", fields = "Version") #"1.2.1"

# Set working environment
setwd("/Users/kibanez/Documents/STRs/ANALYSIS/pipeline_performance/")

# Load golden validation table - EHv3.1.2
#val_data = read.csv("EHv3_avg_VS_EHv3_maxCI/STRVALIDATION_ALLDATA_2019-10-7_ALL_kibanez_EHv312_avg_VS_EHv312_maxCI_checkFXN_withPileup_and_expValidatedData.tsv",
val_data = read.csv("EHv3_avg_VS_EHv3_maxCI/STRVALIDATION_ALLDATA_2019-10-7_ALL_kibanez_EHv312_avg_VS_EHv312_maxCI_checkFXN_251119_withPileup_and_expValidatiedData.tsv",
                    sep = "\t",
                    header = T,
                    stringsAsFactors = F)
dim(val_data)
# 638  24

# Filter the good ones
# 1 - Only keep with `Pileup_quality` == good or Good, `MISSING` and `blanks`
val_data = val_data %>%
  filter(Pileup_quality %in% "Good" | Pileup_quality %in% "good" | Pileup_quality %in% "" | Pileup_quality %in% "MISSING")
dim(val_data)
# 546  24

# 2 - Only keep experimental val numbers (STR_a1, STR_a2) that are integer
val_data = val_data %>%
  filter((!STR_a1 %in% "positive"))
dim(val_data)
# 538  24

val_data = val_data %>%
  filter((!STR_a2 %in% "positive"))
dim(val_data)
# 538  24

val_data = val_data %>%
  filter((!STR_a1 %in% "na"))
dim(val_data)
# 538  24

val_data = val_data %>%
  filter((!STR_a2 %in% "na"))
dim(val_data)
# 521  24

# Let's simplify the data we need from `val_data`
val_data = val_data %>%
  select(LP_Number, locus_bioinfo, STR_a1, STR_a2, EH_a1_avg, EH_a2_avg)

# Just remove all rows having `normal`, `expansion`, `full_mutation`, `premutation`
val_data = val_data %>% filter(!STR_a1 %in% 'normal')
val_data = val_data %>% filter(!STR_a2 %in% 'normal')
val_data = val_data %>% filter(!STR_a1 %in% 'expansion')
val_data = val_data %>% filter(!STR_a2 %in% 'expansion')
val_data = val_data %>% filter(!STR_a1 %in% 'full_mutation')
val_data = val_data %>% filter(!STR_a2 %in% 'full_mutation')
val_data = val_data %>% filter(!STR_a1 %in% 'premutation')
val_data = val_data %>% filter(!STR_a2 %in% 'premutation')

# Remove ONLY COLUMNS (and aligning experimental and EH) in which we do not have a number or we do have `.`
l_dot = which(val_data$STR_a2 == '.')
val_data$EH_a2_avg[l_dot] = 0

l_EXP = which(val_data$STR_a2 == 'EXP')
val_data$EH_a2_avg[l_EXP] = 0

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
  if (val_eh_a2 == 0){
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
# 516  10


# Let's take the important meat: experimentally validated data and EH estimations
exp_alleles_v2 = c(as.integer(val_data$validation_a1), as.integer(val_data$validation_a2))
eh_alleles_v2 = c(as.integer(val_data$eh_a1), as.integer(val_data$eh_a2))
locus_v2 = c(val_data$locus_bioinfo, val_data$locus_bioinfo)

# Remove NAs
index_NA = which(is.na(exp_alleles_v2))
exp_alleles_v2 = exp_alleles_v2[-index_NA]
eh_alleles_v2 = eh_alleles_v2[-index_NA]
locus_v2 = locus_v2[-index_NA]

# Create dataframe with exp, eh, freq for each locus
df_data_with_freq_v2 = data.frame()
l_locus = unique(locus_v2)
for(i in 1:length(l_locus)){
  aux_validation_a1 = val_data %>% filter(locus_bioinfo %in% l_locus[i]) %>% select(validation_a1) %>% pull() %>% as.integer() 
  aux_validation_a2 = val_data %>% filter(locus_bioinfo %in% l_locus[i]) %>% select(validation_a2) %>% pull() %>% as.integer() 
  aux_exp_alleles_v2 = c(aux_validation_a1, aux_validation_a2)
  
  aux_eh_a1 = val_data %>% filter(locus_bioinfo %in% l_locus[i]) %>% select(eh_a1) %>% pull() %>% as.integer() 
  aux_eh_a2 = val_data %>% filter(locus_bioinfo %in% l_locus[i]) %>% select(eh_a2) %>% pull() %>% as.integer() 
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

output_folder = "./figures/EHv3.1.2/"

#df_data_with_freq_v2$exp_alleles[142] = 100

max_value = max(df_data_with_freq_v2$eh_alleles, 
                df_data_with_freq_v2$exp_alleles,
                na.rm = TRUE) + 5

joint_plot = ggplot(df_data_with_freq_v2, 
                    aes(x = eh_alleles, y = exp_alleles)) + 
  geom_point(aes(color = locus, size = number_of_alleles)) + 
  xlim(5,max_value) + 
  ylim(5,max_value) + 
  labs(title = "", 
       x = "Repeat sizes for each allele \n Expansion Hunter (EH-v3.1.2)", 
       y = "Repeat sizes for each allele \n Experimental validation") + 
  geom_abline(method = "lm", formula = y ~ x, linetype = 2, colour = "gray") +  
  coord_equal() +
  guides(size = FALSE)

png("figures/joint_bubble_plot_EHv3.png")
print(joint_plot)
dev.off()

pdf("figures/joint_bubble_plot_EHv3.pdf")
print(joint_plot)
dev.off()

# let's plot each locus independently
for(i in 1:length(l_locus)){
  df_data_with_freq_v2_locus = df_data_with_freq_v2 %>% 
    filter(locus %in% l_locus[i])
  
  max_value = max(df_data_with_freq_v2_locus$eh_alleles, 
                  df_data_with_freq_v2_locus$exp_alleles) + 5
  
  file_name = paste(l_locus[i], "experimental_vs_EHv3.1.2", sep = "_")
  pdf_name = paste(file_name, "pdf", sep = ".")
  png_name = paste(file_name, "png", sep = ".")
  pdf_output = paste(output_folder, pdf_name, sep = "")
  png_output = paste(output_folder, png_name, sep = "")
  
  locus_bubble = ggplot(df_data_with_freq_v2_locus, 
         aes(x = eh_alleles, y = exp_alleles)) + 
    geom_point(aes(color = locus, size = number_of_alleles), show.legend = FALSE) + 
    xlim(5,max_value) + 
    ylim(5,max_value) + 
    labs(title = l_locus[i], 
         x = "Repeat sizes for each allele \n Expansion Hunter (EH-v3.1.2)", 
         y = "Repeat sizes for each allele \n Experimental validation") + 
    theme(plot.title = element_text(hjust = 0.5)) +
    geom_abline(method = "lm", formula = y ~ x, linetype = 2, colour = "gray") +  
    coord_equal()
  
  pdf(pdf_output)
  print(locus_bubble)
  dev.off()
  
  png(png_output)
  print(locus_bubble)
  dev.off()
}

