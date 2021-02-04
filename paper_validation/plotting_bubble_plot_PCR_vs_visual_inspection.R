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
library(cowplot); packageDescription ("cowplot", fields = "Version") #"1.0.0"

# Set working environment
setwd("/Users/kibanez/Documents/STRs/VALIDATION/bubble_plots/")

# Load golden validation table - 509 PCR tests with AT LEAST one allele with exact PCR sizes
val_data = read.csv("./GEL_accuracy_final_not_UCL_considering_PCR_exp_shorter_readLength_270121.tsv",
                    sep = "\t",
                    header = T,
                    stringsAsFactors = F)

dim(val_data)
# 485  10

val_data = read.csv("./GEL_accuracy_final_not_UCL_considering_PCR_exp_larger_readLength_270121.tsv",
                    sep = "\t",
                    header = T,
                    stringsAsFactors = F)

dim(val_data)
# 509  10

output_folder = "./figures/"

# Mike's suggestion - 1
# PCR alleles in X axis
# EHv3 after visual inspection (changing 5 FP and 1 FN values)

# Let's take the important meat
exp_alleles_v2 = c(as.integer(val_data$exp_PCR_a1), as.integer(val_data$exp_PCR_a2))
eh_alleles_v2 = c(as.integer(val_data$EHv312_a1_avg_after_visualQC), as.integer(val_data$EHv312_a2_avg_after_visualQC))
locus_v2 = c(val_data$locus, val_data$locus)

# Remove NAs
index_NA = which(is.na(exp_alleles_v2))
# Remove EXP,NORM,PREMUT,FULL EXP
index_EXP = which(grepl("EXP",exp_alleles_v2))

exp_alleles_v2 = exp_alleles_v2[-index_NA]
eh_alleles_v2 = eh_alleles_v2[-index_NA]
locus_v2 = locus_v2[-index_NA]

# Create dataframe with exp, eh, freq for each locus
df_data_with_freq_v2 = data.frame()
l_locus = unique(locus_v2)
for(i in 1:length(l_locus)){
  #aux_validation_a1 = val_data %>% filter(locus %in% l_locus[i]) %>% select(Truth.Short.Allele) %>% pull() %>% as.integer() 
  #aux_validation_a2 = val_data %>% filter(locus %in% l_locus[i]) %>% select(Truth.Long.Allele) %>% pull() %>% as.integer() 
  aux_validation_a1 = val_data %>% filter(locus %in% l_locus[i]) %>% select(exp_PCR_a1) %>% pull() %>% as.integer() 
  aux_validation_a2 = val_data %>% filter(locus %in% l_locus[i]) %>% select(exp_PCR_a2) %>% pull() %>% as.integer() 
  aux_exp_alleles_v2 = c(aux_validation_a1, aux_validation_a2)
  
  #aux_eh_a1 = val_data %>% filter(locus %in% l_locus[i]) %>% select(min.EHv312.a1) %>% pull() %>% as.integer() 
  #aux_eh_a2 = val_data %>% filter(locus %in% l_locus[i]) %>% select(max.EHv312.a2) %>% pull() %>% as.integer() 
  aux_eh_a1 = val_data %>% filter(locus %in% l_locus[i]) %>% select(EHv312_a1_avg_after_visualQC) %>% pull() %>% as.integer() 
  aux_eh_a2 = val_data %>% filter(locus %in% l_locus[i]) %>% select(EHv312_a2_avg_after_visualQC) %>% pull() %>% as.integer() 
  aux_eh_alleles_v2 = c(aux_eh_a1, aux_eh_a2)
  
  # Since we are keeping alleles for which sometimes an allele has EXP/NORMAL/PREMUT and not the other, we need to remove the same index
  index_na = which(is.na(aux_exp_alleles_v2))
  
  if (length(index_na) > 0){
    data_aux = xyTable(aux_exp_alleles_v2[-index_na], 
                       aux_eh_alleles_v2[-index_na])
  }else{
    data_aux = xyTable(aux_exp_alleles_v2[!is.na(aux_exp_alleles_v2)], 
                       aux_eh_alleles_v2[!is.na(aux_eh_alleles_v2)])
  }
  
  df_data_aux = data.frame(eh_alleles = data_aux$y,
                           exp_alleles = data_aux$x,
                           number_of_alleles = data_aux$number,
                           locus = rep(l_locus[i], length(data_aux$x)))
  
  
  # Concat info per locus
  df_data_with_freq_v2 = rbind(df_data_with_freq_v2,
                               df_data_aux)
  
}

max_value = max(df_data_with_freq_v2$eh_alleles, 
                df_data_with_freq_v2$exp_alleles,
                na.rm = TRUE) + 5
df_strategy1 = df_data_with_freq_v2

# Mike's suggestion - 2
# min and max PCR sizes in X axis
# minEhv3 and maxEHv3 before visual QC in Y axis

# Let's take the important meat
exp_alleles_v2 = c(as.integer(val_data$exp_PCR_a1), as.integer(val_data$exp_PCR_a2))
eh_alleles_v2 = c(as.integer(val_data$EHv312_a1_avg), as.integer(val_data$EHv312_a2_avg))
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
  aux_validation_a1 = val_data %>% filter(locus %in% l_locus[i]) %>% select(exp_PCR_a1) %>% pull() %>% as.integer() 
  aux_validation_a2 = val_data %>% filter(locus %in% l_locus[i]) %>% select(exp_PCR_a2) %>% pull() %>% as.integer() 
  aux_exp_alleles_v2 = c(aux_validation_a1, aux_validation_a2)
  
  aux_eh_a1 = val_data %>% filter(locus %in% l_locus[i]) %>% select(EHv312_a1_avg) %>% pull() %>% as.integer() 
  aux_eh_a2 = val_data %>% filter(locus %in% l_locus[i]) %>% select(EHv312_a2_avg) %>% pull() %>% as.integer() 
  aux_eh_alleles_v2 = c(aux_eh_a1, aux_eh_a2)
  
  index_na = which(is.na(aux_exp_alleles_v2))
  
  if (length(index_na) >0){
    data_aux = xyTable(aux_exp_alleles_v2[-index_na], 
                       aux_eh_alleles_v2[-index_na])
  }else{
    data_aux = xyTable(aux_exp_alleles_v2[!is.na(aux_exp_alleles_v2)], 
                       aux_eh_alleles_v2[!is.na(aux_eh_alleles_v2)])
  }
  
  df_data_aux = data.frame(eh_alleles = data_aux$y,
                           exp_alleles = data_aux$x,
                           number_of_alleles = data_aux$number,
                           locus = rep(l_locus[i], length(data_aux$x)))
  # Concat info per locus
  df_data_with_freq_v2 = rbind(df_data_with_freq_v2,
                               df_data_aux)
  
}

max_value = max(df_data_with_freq_v2$eh_alleles, 
                df_data_with_freq_v2$exp_alleles,
                na.rm = TRUE) + 5

df_strategy2 = df_data_with_freq_v2

group.colors = rainbow(13)

# Combining joint_plot_mike1 and joint_plot_mike2 into a single one
tontz = ggplot() +
geom_point(data = df_strategy2, aes(x = exp_alleles, y = eh_alleles, size = number_of_alleles), color = "#B8B8B8") +
  geom_point(data = df_strategy1, aes(color = factor(locus), x = exp_alleles, y = eh_alleles, size = number_of_alleles), alpha = 0.7) +  
  xlim(5,max_value) +
  ylim(5,max_value) +
  geom_abline(method = "lm", formula = x ~ y, linetype = 2, colour = "gray") +  
  coord_equal() +
  scale_color_manual(values=group.colors) +
  labs(title = "", 
       y = "EH repeat sizes", 
       x = "PCR repeat sizes") + 
  theme_light() +
  theme(legend.title = element_blank(),
        text = element_text(size=13),
        axis.text.x.top = element_text()) +
  guides(size = FALSE) 

png("./figures/FigureS3_418PCRtests_filtering_out_NCL_shorterThanReadLengthLANCET_600dpi_040221.png",units="in", width=5, height=5, res=600)
png("./figures/FigureS3_418PCRtests_filtering_out_NCL_largerThanReadLengthLANCET_600dpi_040221.png",units="in", width=5, height=5, res=600)
print(tontz)
dev.off()

# breakdown by locus
# First enrich `val_data` with premutation cut-offs
df_cutoffs = data.frame(locus = c("AR", "ATN1", "ATXN1", "ATXN2", "ATXN3", "ATXN7", "CACNA1A", "C9orf72", "DMPK", "HTT", "FMR1", "FXN", "TBP"),
                        premut_cutoff = c(34,34,35,31,43,34,17,30,50,35,55,44,41),
                        stringsAsFactors = F)

val_data = left_join(val_data,
                     df_cutoffs,
                     by = "locus")
df_strategy1 = left_join(df_strategy1,
                         df_cutoffs,
                         by = "locus")
df_strategy2 = left_join(df_strategy2,
                         df_cutoffs,
                         by = "locus")

df_strategy1$locus = as.factor(df_strategy1$locus)
df_strategy2$locus = as.factor(df_strategy2$locus)

all_except_DMPK = c("AR", "ATN1", "ATXN1", "ATXN2", "ATXN3", "ATXN7", "CACNA1A", "C9orf72", "HTT", "FMR1", "FXN", "TBP")
breakdown_by_locus = ggplot(df_strategy1) +
  geom_point(data = df_strategy2, aes(x = exp_alleles, y = eh_alleles, size = number_of_alleles), color = "#B8B8B8") +
  geom_point(data = df_strategy1, aes(color = factor(locus), x = exp_alleles, y = eh_alleles, size = number_of_alleles), alpha = 0.7) +  
  geom_vline(data = df_strategy1, aes(group = locus, xintercept=as.numeric(premut_cutoff)), color ="red", lwd=0.3, lty=4) +
  geom_hline(data = df_strategy1, aes(group = locus, yintercept=as.numeric(premut_cutoff)), color ="red", lwd=0.3, lty=4) +
  
  xlim(5,max_value) +
  ylim(5,max_value) +
  geom_abline(method = "lm", formula = x ~ y, linetype = 2, colour = "gray") +  
  coord_equal() +
  labs(title = "", 
       y = "EH repeat sizes", 
       x = "PCR repeat sizes") + 
  scale_color_manual(values=group.colors) +
  theme(legend.title = element_blank(),
        axis.text.x.top = element_text()) +
  guides(size = FALSE, color = FALSE) +
  facet_wrap(locus~ .) 
 
  

#png("./figures/Figure2B_LANCET_filter_all_NCL_shorterThanReadLength_600dpi_190121_withDMPK.png",units="in", width=5, height=5, res=600)
png("./figures/Figure2B_LANCET_filter_all_NCL_largerThanReadLength_600dpi_190121.png",units="in", width=5, height=5, res=600)
print(breakdown_by_locus)
dev.off()


# For repeats larger than the read-length, considering all repeat-sizes
# loading the 2nd `val_data`
# Let's split in two, before 150 and after
# load colour per locus
all_loci = c("AR", "ATN1", "ATXN1", "ATXN2", "ATXN3", "ATXN7", "C9orf72", "CACNA1A", "DMPK", "FMR1", "FXN", "HTT", "TBP")

all_loci = c("AR", "ATN1", "ATXN1", "ATXN2", "ATXN3", "ATXN7", "C9orf72", "CACNA1A", "DMPK - repeats ≤ read-length", "FMR1", "FXN - repeats ≤ read-length", "HTT", "TBP")
df_strategy1$locus = gsub("DMPK - repeats > read-length", "DMPK - repeats ≤ read-length", df_strategy1$locus)
df_strategy2$locus = gsub("DMPK - repeats > read-length", "DMPK - repeats ≤ read-length", df_strategy2$locus)
df_strategy1$locus = gsub("FXN - repeats > read-length", "FXN - repeats ≤ read-length", df_strategy1$locus)
df_strategy2$locus = gsub("FXN - repeats > read-length", "FXN - repeats ≤ read-length", df_strategy2$locus)

for (locus_name in 1:length(all_loci)){
  colour_locus = group.colors[locus_name]
  
  # Define first max_value
  max_value_df1 = max(df_strategy1 %>% filter(locus %in% all_loci[locus_name], exp_alleles <=150) %>% select(eh_alleles) %>% pull(), 
                      df_strategy1 %>% filter(locus %in% all_loci[locus_name], exp_alleles <=150) %>% select(exp_alleles) %>% pull(),
                  na.rm = TRUE) + 5
  
  max_value_df2 = max(df_strategy2 %>% filter(locus %in% all_loci[locus_name], exp_alleles <=150) %>% select(eh_alleles) %>% pull(), 
                      df_strategy2 %>% filter(locus %in% all_loci[locus_name], exp_alleles <=150) %>% select(exp_alleles) %>% pull(),
                      na.rm = TRUE) + 5
  
  max_value = max(max_value_df1,
                  max_value_df2)
  
  breakdown_by_locus_shorter = ggplot(df_strategy1 %>% filter(locus %in% all_loci[locus_name], exp_alleles <= 150)) +
    geom_point(data = df_strategy2 %>% filter(locus %in% all_loci[locus_name], exp_alleles <=150), aes(x = exp_alleles, y = eh_alleles, size = number_of_alleles), color = "#B8B8B8") +
    geom_point(data = df_strategy1 %>% filter(locus %in% all_loci[locus_name], exp_alleles <=150), aes(x = exp_alleles, y = eh_alleles, size = number_of_alleles), color = colour_locus, alpha = 0.7) +  
    geom_vline(data = df_strategy1 %>% filter(locus %in% all_loci[locus_name], exp_alleles <=150), aes(group = locus, xintercept=as.numeric(premut_cutoff)), color ="red", lwd=0.3, lty=4) +
    geom_hline(data = df_strategy1 %>% filter(locus %in% all_loci[locus_name], exp_alleles <=150), aes(group = locus, yintercept=as.numeric(premut_cutoff)), color ="red", lwd=0.3, lty=4) +
    xlim(5,max_value) +
    ylim(5,max_value) +
    geom_abline(method = "lm", formula = x ~ y, linetype = 2, colour = "gray") +  
    coord_equal() +
    labs(title = "", 
         y = "", 
         x = "") + 
#    scale_color_manual(values=colour_locus) +
    theme(legend.title = element_blank(),
          #axis.text.x.top = element_text(),
          axis.title.x = element_blank(),
          axis.text.x = element_text(size=25, angle = 45, vjust=1, hjust = 1),
          axis.text.y = element_text(size=25, angle = 45, vjust=0, hjust = 1)) +
          #text=element_text(size=35),
          #strip.text = element_text(size=18)) +
    guides(size = FALSE, color = FALSE) 
    #facet_wrap(locus~ .) 
  
  # Define first max_value
  max_value_df1 = max(df_strategy1 %>% filter(locus %in% all_loci[locus_name]) %>% select(eh_alleles) %>% pull(), 
                      df_strategy1 %>% filter(locus %in% all_loci[locus_name]) %>% select(exp_alleles) %>% pull(),
                      na.rm = TRUE) + 5
  
  max_value_df2 = max(df_strategy2 %>% filter(locus %in% all_loci[locus_name]) %>% select(eh_alleles) %>% pull(), 
                      df_strategy2 %>% filter(locus %in% all_loci[locus_name]) %>% select(exp_alleles) %>% pull(),
                      na.rm = TRUE) + 5
  
  max_value = max(max_value_df1,
                  max_value_df2)
  
  breakdown_by_locus_shorter_merged = ggplot(df_strategy1 %>% filter(locus %in% all_loci[locus_name])) +
    geom_point(data = df_strategy2 %>% filter(locus %in% all_loci[locus_name]), aes(x = exp_alleles, y = eh_alleles, size = number_of_alleles), color = "#B8B8B8") +
    geom_point(data = df_strategy1 %>% filter(locus %in% all_loci[locus_name]), aes(x = exp_alleles, y = eh_alleles, size = number_of_alleles), color = colour_locus, alpha = 0.7) +  
    geom_vline(data = df_strategy1 %>% filter(locus %in% all_loci[locus_name]), aes(group = locus, xintercept=as.numeric(premut_cutoff)), color ="red", lwd=0.3, lty=4) +
    geom_hline(data = df_strategy1 %>% filter(locus %in% all_loci[locus_name]), aes(group = locus, yintercept=as.numeric(premut_cutoff)), color ="red", lwd=0.3, lty=4) +
    
    xlim(5,max_value) +
    ylim(5,max_value) +
    geom_abline(method = "lm", formula = x ~ y, linetype = 2, colour = "gray") +  
    coord_equal() +
    labs(title = "", 
         y = "", 
         x = "") + 
    #scale_color_manual(values=group.colors) +
    theme(legend.title = element_blank(),
          #axis.text.x.top = element_text(),
          axis.title.x = element_blank(),
          axis.text.x = element_text(size=25, angle = 45, vjust=1, hjust = 1),
          axis.text.y = element_text(size=25, angle = 45, vjust=0, hjust = 1)) +
    guides(size = FALSE, color = FALSE) 
  
  breakdown_by_locus_larger = ggplot(df_strategy1 %>% filter(locus %in% all_loci[locus_name], exp_alleles > 150)) +
    geom_point(data = df_strategy2 %>% filter(locus %in% all_loci[locus_name], exp_alleles >150), aes(x = exp_alleles, y = eh_alleles, size = number_of_alleles), color = "#B8B8B8") +
    geom_point(data = df_strategy1 %>% filter(locus %in% all_loci[locus_name], exp_alleles >150), aes(x = exp_alleles, y = eh_alleles, size = number_of_alleles), color = colour_locus, alpha = 0.7) +  
    geom_vline(data = df_strategy1 %>% filter(locus %in% all_loci[locus_name], exp_alleles >150), aes(group = locus, xintercept=as.numeric(premut_cutoff)), color ="red", lwd=0.3, lty=4) +
    geom_hline(data = df_strategy1 %>% filter(locus %in% all_loci[locus_name], exp_alleles >150), aes(group = locus, yintercept=as.numeric(premut_cutoff)), color ="red", lwd=0.3, lty=4) +
    
    xlim(5,max_value) +
    ylim(5,max_value) +
    geom_abline(method = "lm", formula = x ~ y, linetype = 2, colour = "gray") +  
    coord_equal() +
    labs(title = "", 
         y = "", 
         x = "") + 
    #scale_color_manual(values=group.colors) +  
    theme(legend.title = element_blank(),
          #axis.text.x.top = element_text(),
          axis.title.x = element_blank(),
          axis.text.x = element_text(size=25, angle = 45, vjust=1, hjust = 1),
          axis.text.y = element_text(size=25, angle = 45, vjust=0, hjust = 1)) +
    guides(size = FALSE, color = FALSE) 
  
  file_name = all_loci[locus_name]
  file_name = paste("./figures/", file_name, sep = "")
  file_name_short = paste(file_name, "shorterThanReadLength_600dpi_040221.png", sep = "_") 
  file_name_large = paste(file_name, "largerThanReadLength_600dpi_040221.png", sep = "_") 
  file_name_merged = paste(file_name, "merged_600dpi_04020121.png", sep = "_") 
  
  png(file_name_short,units="in", width=5, height=5, res=600)
  print(breakdown_by_locus_shorter)
  dev.off()
  
  png(file_name_merged,units="in", width=5, height=5, res=600)
  print(breakdown_by_locus_shorter_merged)
  dev.off()
  
  png(file_name_large,units="in", width=9, height=6, res=600)
  print(breakdown_by_locus_larger)
  dev.off()
}

