# Objective: explore the distribution of STR allele repeat-sizes across the validation golden table dataset
date ()
Sys.info ()[c("nodename", "user")]
commandArgs ()
rm (list = ls ())
R.version.string ## "R version 3.5.3 (2019-03-11)

library(ggplot2); packageDescription ("ggplot2", fields = "Version") #"2.2.1"
library(grid); packageDescription ("grid", fields = "Version") #"3.3.2"
library(gridExtra); packageDescription ("gridExtra", fields = "Version") #"2.2.1"
library(reshape2); packageDescription ("reshape2", fields = "Version") #"1.4.2"
require(dplyr); packageDescription ("dplyr", fields = "Version") #"0.7.1"

# Set working directory
setwd("~/Documents/STRs/ANALYSIS/EHdn/")

# EHdn
# We want to see the distribution of the number of repeats / WGS estimated by EHdn
l_ehdn_b37 = list.files("~/Documents/STRs/ANALYSIS/EHdn/output_EHdn_v0.8.0/numberRepeats/GRCh37/")
l_ehdn_b38 = list.files("~/Documents/STRs/ANALYSIS/EHdn/output_EHdn_v0.8.0/numberRepeats/GRCh38/")

l_ehdn_b37 = paste("~/Documents/STRs/ANALYSIS/EHdn/output_EHdn_v0.8.0/numberRepeats/GRCh37/", l_ehdn_b37, sep = "")
l_ehdn_b38 = paste("~/Documents/STRs/ANALYSIS/EHdn/output_EHdn_v0.8.0/numberRepeats/GRCh38/", l_ehdn_b38, sep = "")


## GRCh37
df_eh_b37 = data.frame()
for (i in 1:length(l_ehdn_b37)){
  file_name = strsplit(l_ehdn_b37[i], "/")[[1]][9]
  id_name = strsplit(file_name, "\\.")[[1]][1]
  number_repeats = length(readLines(l_ehdn_b37[i]))
  df_eh_b37 = rbind(df_eh_b37, cbind(id_name, number_repeats))
}
dim(df_eh_b37)
# 87  2

df_eh_b37$id_name = as.character(df_eh_b37$id_name)
df_eh_b37$number_repeats = as.numeric(as.character(df_eh_b37$number_repeats))

png("./figures/plot_numberRepeats_EHdn-v0.8.0_WGS_GRCh37.png")
ggplot(unique(df_eh_b37), aes(x = id_name, y = number_repeats)) + 
  geom_bar(stat = "identity") + 
  ylab("number of repeats estimated by EHdn-v0.8.0") + 
  xlab("87 validation WGS samples - GRCh37") 
dev.off()


## GRCh38
df_eh_b38 = data.frame()
for (i in 1:length(l_ehdn_b38)){
  file_name = strsplit(l_ehdn_b38[i], "/")[[1]][9]
  id_name = strsplit(file_name, "\\.")[[1]][1]
  number_repeats = length(readLines(l_ehdn_b38[i]))
  df_eh_b38 = rbind(df_eh_b38, cbind(id_name, number_repeats))
}
dim(df_eh_b38)
# 127  2


df_eh_b38$id_name = as.character(df_eh_b38$id_name)
df_eh_b38$number_repeats = as.numeric(as.character(df_eh_b38$number_repeats))

png("figures/plot_numberRepeats_EHdn-v0.8.0_WGS_GRCh38.png")
ggplot(unique(df_eh_b38), aes(x = id_name, y = number_repeats)) + 
  geom_bar(stat = "identity") + 
  ylab("number of repeats estimated by EHdn-v0.8.0") + 
  xlab("127 validation WGS samples - GRCh38") 
dev.off()


# Let's plot GRCh37 and GRCh38 together and see whether there is a difference on EHdn performance depending on the genome build
df_eh_b37$genome = rep("GRCh37", length(df_eh_b37$id_name))
df_eh_b38$genome = rep("GRCh38", length(df_eh_b38$id_name))

df_ehdn = rbind(df_eh_b37,
                df_eh_b38)
dim(df_ehdn)
# 214  3

png("figures/plot_numberRepeats_EHdn-v0.8.0_WGS_GRCh37_and_GRCh38.png")
ggplot(df_ehdn, aes(x = id_name, y = number_repeats, color = genome)) +
  geom_bar(stat = "identity") + 
  ylab("number of repeats estimated by EHdn-v0.8.0") + 
  xlab("214 validation WGS samples") 
dev.off()
