# Objective: compare GEL validation golden table comparing EHv2.5.5 vs EHv3.1.2
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
setwd("/Users/kibanez/Documents/STRs/ANALYSIS/pipeline_performance/GEL_val_data_EHv255_vs_EHv312/")

# Load golden validation table - EHv2.5.5
val_data_v2 = read.csv("../EHv2_avg_VS_EHv2_maxCI/STRVALIDATION_ALLDATA_2019-10-7_ALL_kibanez_EHv255_avg_VS_EHv255_maxCI_checkFXN_withPileup_and_expValidatedData_tweaking_ATN1_updated_AR_from_NHNN.txt",
                    sep = "\t",
                    header = T,
                    stringsAsFactors = F)
dim(val_data_v2)
# 638  20

# Load golden validation table - EHv3.1.2
val_data_v3 = read.csv("../EHv3_avg_VS_EHv3_maxCI/STRVALIDATION_ALLDATA_2019-10-7_ALL_kibanez_EHv312_avg_VS_EHv312_maxCI_ClassiByAllele.tsv",
                    sep = "\t",
                    header = T,
                    stringsAsFactors = F)
dim(val_data_v3)
# 638  19

# Merge V2 and V3, after creating a new column with the EHversion
val_data_v2$EH_version = rep("EH-v2.5.5", length(val_data_v2$locus_bioinfo))
val_data_v3$EH_version = rep("EH-v3.1.2", length(val_data_v3$locus_bioinfo))

# Merge
col_intersected = intersect(colnames(val_data_v2), colnames(val_data_v3))
val_data_v2 = val_data_v2 %>% select(col_intersected)
val_data_v3 = val_data_v3 %>% select(col_intersected)

val_data = rbind(val_data_v2,
                 val_data_v3)
dim(val_data)
# 1276  16

# Let's simplify the data we need from `val_data`
val_data = val_data %>%
  select(LP_Number, locus_bioinfo, EH_a1_avg, EH_a2_avg, EH_version)

# Before creating plots, estimate the frequency of each allele repeat-size
# Create dataframe with exp, eh, freq for each locus
df_data_with_freq = data.frame()
l_locus = unique(val_data$locus_bioinfo)
for(i in 1:length(l_locus)){
  aux_ehv2_a1 = val_data %>% filter(locus_bioinfo %in% l_locus[i], EH_version %in% "EH-v2.5.5") %>% select(EH_a1_avg) %>% pull() %>% as.integer() 
  aux_ehv2_a2 = val_data %>% filter(locus_bioinfo %in% l_locus[i], EH_version %in% "EH-v2.5.5") %>% select(EH_a2_avg) %>% pull() %>% as.integer() 
  aux_ehv2_alleles = c(aux_ehv2_a1, aux_ehv2_a2)
  
  aux_ehv3_a1 = val_data %>% filter(locus_bioinfo %in% l_locus[i], EH_version %in% "EH-v3.1.2") %>% select(EH_a1_avg) %>% pull() %>% as.integer() 
  aux_ehv3_a2 = val_data %>% filter(locus_bioinfo %in% l_locus[i], EH_version %in% "EH-v3.1.2") %>% select(EH_a2_avg) %>% pull() %>% as.integer() 
  aux_ehv3_alleles = c(aux_ehv3_a1, aux_ehv3_a2)
  
  data_aux = xyTable(aux_ehv2_alleles, aux_ehv3_alleles)
  
  df_data_aux = data.frame(ehv3_alleles = data_aux$y,
                           ehv2_alleles = data_aux$x,
                           number_of_alleles = data_aux$number,
                           locus = rep(l_locus[i], length(data_aux$x)))
  # Concat info per locus
  df_data_with_freq = rbind(df_data_with_freq,
                               df_data_aux)
}

output_folder = "./figures/"

max_value = max(df_data_with_freq$ehv2_alleles, 
                df_data_with_freq$ehv3_alleles,
                na.rm = TRUE) + 5



# Joint all loci together
joint_plot = ggplot(df_data_with_freq, 
                    aes(x = ehv3_alleles, y = ehv2_alleles)) + 
  geom_point(aes(color = locus, size = number_of_alleles)) + 
  xlim(5,max_value) + 
  ylim(5,max_value) + 
  labs(title = "", 
       y = "Repeat sizes for each allele \n Expansion Hunter (EH-v3.1.2)", 
       x = "Repeat sizes for each allele \n Expansion Hunter (EH-v2.5.5)") + 
  geom_abline(method = "lm", formula = x ~ y, linetype = 2, colour = "gray") +  
  coord_equal() +
  guides(size = FALSE)

png("figures/joint_bubble_plot_GEL_golden_val_table_EHv2_VS_EHv3.png", units="in", width=5, height=5, res=300)
print(joint_plot)
dev.off()

pdf("figures/joint_bubble_plot_GEL_golden_val_table_EHv2_VS_EHv3.pdf")
print(joint_plot)
dev.off()




