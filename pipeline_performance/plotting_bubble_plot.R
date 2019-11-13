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

# Load golden validation table - EHv2.5.5
val_data = read.csv("EHv2_avg_VS_EHv2_maxCI/STRVALIDATION_ALLDATA_2019-10-7_ALL_kibanez_EHv255_avg_VS_EHv255_maxCI_checkFXN_withPileup_and_expValidatedData.tsv",
                    sep = "\t",
                    header = T,
                    stringsAsFactors = F)
dim(val_data)
# 660  19

# Filter the good ones
# 1 - Only keep with `Pileup_quality` == good or Good
val_data = val_data %>%
  filter(Pileup_quality %in% "Good" | Pileup_quality %in% "good")
dim(val_data)
# 443  19

# 2 - Only keep experimental val numbers (STR_a1, STR_a2) that are integer
val_data = val_data %>%
  filter((!STR_a1 %in% "positive"))
dim(val_data)
# 436  19 

val_data = val_data %>%
  filter((!STR_a2 %in% "positive"))
dim(val_data)
# 436  19 

val_data = val_data %>%
  filter((!STR_a1 %in% "na"))
dim(val_data)
# 436  19 

val_data = val_data %>%
  filter((!STR_a2 %in% "na"))
dim(val_data)
# 419  19 




exp_alleles = as.integer(df_locus_final_wessex$validation)
eh_alleles = as.integer(df_locus_final_wessex$eh)

data_with_freq = xyTable(exp_alleles, eh_alleles)
data_with_freq2 = data.frame(eh_alleles = data_with_freq$y, exp_alleles = data_with_freq$x, number_of_alleles = data_with_freq$number)

max_value = max(data_with_freq2$eh_alleles, data_with_freq2$exp_alleles) + 5
png(paste("Correlation_EH_experimental_merging_delivery_versions_WESSEX", paste(l_loci[i],"png", sep = '.'), sep = "_"))
ggplot(data_with_freq2, aes(x = eh_alleles, y = exp_alleles, size = number_of_alleles)) + geom_point(alpha = 0.7, colour = "blue") + xlim(5,max_value) + ylim(5,max_value) + labs(title = paste("Correlation on repeat sizes: EH vs experimental validation", l_loci[i], sep=' '), x = "Repeat sizes for each allele \n Expansion Hunter", y = "Repeat sizes for each allele \n Experimental validation") + geom_abline(method = "lm", formula = y ~ x, linetype = 2, colour = "gray") +  coord_equal()
dev.off()  

