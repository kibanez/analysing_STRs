# Objective: bubble plot for HTT (poster)
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
library(RColorBrewer); packageDescription ("RColorBrewer", fields = "Version") #"1.1-2"

# Set working environment
setwd("/Users/kibanez/Documents/STRs/VALIDATION/")

# Load golden validation table - EHv2.5.5
val_data = read.csv("HTT_validation.tsv",
                    sep = "\t",
                    header = T,
                    stringsAsFactors = F)
dim(val_data)
# 76  6

# Define a1 as min value of both alleles: for PCR and EH
val_data = val_data %>%
  group_by(LP_number) %>%
  mutate(new_PCR_a1 = min(exp_PCR_a1, exp_PCR_a2),
         new_PCR_a2 = max(exp_PCR_a1, exp_PCR_a2),
         new_EH_a1 = min(EHv255_a1_avg, EHv255_a2_avg),
         new_EH_a2 = max(EHv255_a1_avg, EHv255_a2_avg)) %>%
  ungroup() %>%
  as.data.frame()

val_data = val_data %>%
  select(locus, LP_number, new_PCR_a1, new_PCR_a2, new_EH_a1, new_EH_a2)

output_folder = "./figures/"

exp_alleles = c(val_data$new_PCR_a1, val_data$new_PCR_a2)
eh_alleles = c(val_data$new_EH_a1, val_data$new_EH_a2)

dena = data.frame(exp_data = exp_alleles,
                  eh_alleles = eh_alleles)

data_aux = xyTable(exp_alleles, eh_alleles)

df_data_aux = data.frame(eh_alleles = data_aux$y,
                         exp_alleles = data_aux$x,
                         number_of_alleles = data_aux$number)

max_value = max(df_data_aux$eh_alleles, df_data_aux$exp_alleles)                         

htt_bubble = ggplot(df_data_aux, 
                    aes(x = exp_alleles, y = eh_alleles, colour = "red")) + 
  geom_point(aes(fill = "red", size = number_of_alleles)) + 
  xlim(5,max_value) + 
  ylim(5,max_value) + 
  labs(title = "", 
       y = "Repeat sizes for each allele - smaller or eq to 70 repeats \n Expansion Hunter (EH-v2.5.5)", 
       x = "Repeat sizes for each allele \n Experimental validation") + 
  geom_abline(method = "lm", formula = x ~ y, linetype = 2, colour = "gray") +  
  coord_equal() +
  theme(legend.title = element_blank()) + 
  guides(size = FALSE)


png("./figures/HTT_bubble_plot_PCR_vs_EHv2.5.5.png", units="in", width=5, height=5, res=300)
print(htt_bubble)
dev.off()
