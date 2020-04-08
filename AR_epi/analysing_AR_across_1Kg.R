# Objective: analyse AR epi across EUR genomes in the 1Kg project
date ()
Sys.info ()[c("nodename", "user")]
commandArgs ()
rm (list = ls ())
R.version.string ## "R version 3.6.1 (2019-07-05)"

library(ggplot2); packageDescription ("ggplot2", fields = "Version") #"3.2.1"
library(grid); packageDescription ("grid", fields = "Version") #"3.6.2"
library(gridExtra); packageDescription ("gridExtra", fields = "Version") #"2.3"
library(reshape2); packageDescription ("reshape2", fields = "Version") #"1.4.3"
require(dplyr); packageDescription ("dplyr", fields = "Version") #"0.8.3"

# Set working dir
setwd("~/Documents/STRs/ANALYSIS/AR_kennedy/1Kg/")

# Load 1Kg population index data
popu_info = read.csv("~/Documents/STRs/ANALYSIS/population_research/1kg/integrated_call_samples_v2.20130502.ALL.ped",
                     sep = "\t",
                     header = T,
                     stringsAsFactors = F)
dim(popu_info)
# 3691  17

# Load EHv3.2.2 STR merged data for each sub-population
df_merged = data.frame()
l_popus = unique(popu_info$Population)

# Remove CHD from `l_popus`
l_popus = l_popus[-19]

# Define super-populations
l_superpopu = c("AFR", "AMR", "EAS", "EUR", "SAS")

# Define sub-population and super-population
superpopulations = c("AFR","AFR","AFR","AFR","AFR","AFR","AFR", 
                     "AMR", "AMR","AMR","AMR",
                     "EUR","EUR","EUR","EUR","EUR", 
                     "EAS","EAS","EAS","EAS","EAS",
                     "SAS","SAS","SAS","SAS","SAS")
sub_populations = c("ESN", "YRI", "GWD", "LWK", "MSL", "ACB", "ASW", 
                    "MXL", "PUR", "PEL", "CLM",
                    "IBS", "TSI", "GBR", "CEU", "FIN",
                    "JPT", "CHS", "CHB", "CDX", "KHV",
                    "PJL", "STU", "BEB", "ITU", "GIH")

popu_1kg = data.frame(cbind(superpopulations, sub_populations))
popu_1kg$superpopulations = as.character(popu_1kg$superpopulations)
popu_1kg$sub_populations = as.character(popu_1kg$sub_populations)

for (i in 1:length(l_popus)){
  popu_aux = paste("~/Documents/STRs/ANALYSIS/population_research/1Kg/data/", l_popus[i] ,sep = "")
  file_aux = list.files(paste(popu_aux, "merged", sep = "/"))
  file_aux = paste(paste(popu_aux, "merged", sep = "/"), file_aux, sep = "/")
  
  df_aux = read.csv(file_aux,
                    sep  = "\t",
                    stringsAsFactors = F,
                    header = T)
  
  index_subpopu = match(l_popus[i], popu_1kg$sub_populations)
  superpopu = popu_1kg$superpopulations[index_subpopu]
  
  df_aux = df_aux %>% 
    mutate(population = l_popus[i],
           superpopulation = superpopu)
  
  df_merged = rbind(df_merged,
                    df_aux)
  
}

dim(df_merged)
# 14107  14

# Let's select `AR` locus
ar_merged = df_merged %>%
  filter(gene %in% "AR") 
ar_merged = unique(ar_merged)
dim(ar_merged)
# 413  14

# Interesting thing: there are 4 genomes with an expanded AR expansion
ar_merged %>% filter(allele > 35)
#gene allele num_samples
#1   AR     36           1
#2   AR     39           1
#3   AR     38           1
#4   AR     37           1
 


# Plot percentage of each allele within each subpopulation
violin_plot = ggplot(eur_merged, aes(x = subpopu, y=allele, fill = subpopu)) +
  geom_violin() +
  coord_flip() +
  xlab("Sub-population across EUR in 1Kg dataset") + 
  ylab("Repeat sizes (repeat units)") + 
  ggtitle("AR across EUR sub-populations within 1Kg") +
  geom_boxplot(width=0.1) +
  guides(fill=guide_legend(reverse=TRUE)) 

violin_output_png = "./plots/violin_plots_across_EUR_subpopus.png"
png(violin_output_png, units="in", width=5, height=5, res=300)
print(violin_plot)
dev.off()

# Plot joint ancestry repeat-size distributions
min_value = min(eur_merged$allele, 34)
max_value = max(eur_merged$allele, 38)
joint_plot = ggplot(eur_merged, aes(x = allele, y = percent_subpopu, group = subpopu, color = subpopu)) + 
  geom_line() + 
  ylab("Allele frequency") + 
  xlab("Repeat sizes (repeat units)") + 
  geom_vline(xintercept = 34, colour = 'blue', lty = 2) + 
  geom_vline(xintercept = 38, colour = 'red', lty = 2) + 
  coord_cartesian(xlim = c(min_value,max_value))

joint_ancestries = "./plots/joint_subpopus_in_AR_1Kg.png"
png(joint_ancestries, units="in", width=5, height=5, res=300)
print(joint_plot)
dev.off()
