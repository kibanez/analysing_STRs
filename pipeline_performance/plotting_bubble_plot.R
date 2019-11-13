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

# Let's simplify the data we need from `val_data`
val_data = val_data %>%
  select(LP_Number, locus_bioinfo, STR_a1, STR_a2, EH_a1_avg, EH_a2_avg)

# in case there is only 1 allele, we should have 1 allele for expValidation and EH
# in case there are 2 alleles, then the `a1` for expValidation and EH should be the minimum (minor repeat size)
for (i in 1:length(val_data$LP_Number)){
  val_validation_a1 = val_data$STR_a1[i]
  val_validation_a2 = val_data$STR_a2[i]
  val_eh_a1 = val_data$EH_a1_avg[i]
  val_eh_a2 = val_data$EH_a2_avg[i]
  
  # For alleles that there is no estimation (EH only calls 1 allele) we do have a `0` -- we ignore/avoid this step in these cases
  if (val_eh_a2 == "0"){
    min_validation = val_validation_a1
    min_eh = val_eh_a1
    max_validation = '.'
    max_eh = '.'
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
# 419  10

# Exclude `validation_a2` or `eh_a2` == '.'
val_data = val_data %>%
  filter(!validation_a2 %in% '.')
val_data = val_data %>%
  filter(!eh_a2 %in% '.')
dim(val_data)
# 417  10

# Let's take the important meat: experimentally validated data and EH estimations
exp_alleles_v2 = c(as.integer(val_data$STR_a1), as.integer(val_data$STR_a2))
eh_alleles_v2 = c(as.integer(val_data$EH_a1_avg), as.integer(val_data$EH_a2_avg))

data_with_freq_v2 = xyTable(exp_alleles_v2, eh_alleles_v2)
df_data_with_freq_v2 = data.frame(eh_alleles = data_with_freq_v2$y, 
                                  exp_alleles = data_with_freq_v2$x, 
                                  number_of_alleles = data_with_freq_v2$number)

max_value = max(df_data_with_freq_v2$eh_alleles, 
                df_data_with_freq_v2$exp_alleles) + 5

png(paste("Correlation_EH_experimental_merging_delivery_versions_WESSEX", paste(l_loci[i],"png", sep = '.'), sep = "_"))
ggplot(df_data_with_freq_v2, 
       aes(x = eh_alleles, y = exp_alleles, size = number_of_alleles)) + 
  geom_point(alpha = 0.7, colour = "blue") + 
  xlim(5,max_value) + 
  ylim(5,max_value) + 
  labs(title = paste("Correlation on repeat sizes: EH vs experimental validation", l_loci[i], sep=' '), x = "Repeat sizes for each allele \n Expansion Hunter", y = "Repeat sizes for each allele \n Experimental validation") + 
  geom_abline(method = "lm", formula = y ~ x, linetype = 2, colour = "gray") +  
  coord_equal()
dev.off()  

