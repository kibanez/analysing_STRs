# Objective: analyse loci shared by Andrew Sharp (i.e. SHARP genes) within 100kGP - batch august 2020
# For EHv322
date ()
Sys.info ()[c("nodename", "user")]
commandArgs ()
rm (list = ls ())
R.version.string ## "R version 3.6.1 (2019-07-05)"

# libraries
library(ggplot2); packageDescription ("ggplot2", fields = "Version") #"3.2.1"
library(grid); packageDescription ("grid", fields = "Version") #"3.6.1"
library(gridExtra); packageDescription ("gridExtra", fields = "Version") #"2.3"
library(reshape2); packageDescription ("reshape2", fields = "Version") #"1.4.3"
require(dplyr); packageDescription ("dplyr", fields = "Version") #"0.8.3"
library(magick); packageDescription ("magick", fields = "Version") #"2.4.0"
library(cowplot); packageDescription ("cowplot", fields = "Version") #"1.0.0"
library(ggpubr); packageDescription ("ggpubr", fields = "Version") #"1.0.0"

# Set working dir
setwd("~/Documents/STRs/ANALYSIS/SHARP/")

# load merged august data
merged_data = read.csv("~/Documents/STRs/data/research/batch_august2020/output_EHv3.2.2_vcfs/merged/merged_93446_genomes_EHv322_batch_august2020.tsv",
                       stringsAsFactors = F,
                       sep = "\t",
                       header = T)
dim(merged_data)
# 27238  12

# load clinical data
clin_data = read.csv("~/Documents/STRs/clinical_data/clinical_research_cohort/clinical_data_research_cohort_93614_PIDs_merging_RE_V1toV10.tsv",
                     sep = "\t",
                     stringsAsFactors = F,
                     header = T)
dim(clin_data)
# 152337  32

# List of platekeys corresponding to ONLY PROBANDS
df_only_probands = clin_data %>%
  filter(is.na(biological_relationship_to_proband) |
           biological_relationship_to_proband %in% "N/A" | 
           biological_relationship_to_proband %in% "Proband" |
           programme %in% "Cancer")

l_platekeys_probands = df_only_probands %>%
  select(list_platekeys1) %>%
  unique() %>%
  pull()
length(l_platekeys_probands)
# 52191

# There are some platekeys (16k) that have ',', which means that PID is associated with more than one platekey
l_platekeys_probands_unique = c()
for (i in 1:length(l_platekeys_probands)){
  if (grepl(',',l_platekeys_probands[i])){
    list_platekeys = strsplit(l_platekeys_probands[i], ",")[[1]]
    list_platekeys = gsub(" ", "", list_platekeys, fixed = TRUE)
    l_platekeys_probands_unique = c(l_platekeys_probands_unique,
                                    max(list_platekeys))
  }else{
    l_platekeys_probands_unique = c(l_platekeys_probands_unique,
                                    l_platekeys_probands[i])
  }
}
length(l_platekeys_probands_unique)
# 52191

# List of platekeys corresponding to ONLY PROBANDS but NOT in Neuro
# First probands
df_only_probands_notNeuro = df_only_probands %>%
  filter(!grepl("neuro", list_disease_group, ignore.case = TRUE))
dim(df_only_probands_notNeuro)
# 63266 32

l_platekeys_probands_notNeuro = df_only_probands_notNeuro %>%
  select(list_platekeys1) %>%
  unique() %>%
  pull()
length(l_platekeys_probands_notNeuro)
# 37701

# There are some platekeys (16k) that have ',', which means that PID is associated with more than one platekey
l_platekeys_probands_notNeuro_unique = c()
for (i in 1:length(l_platekeys_probands_notNeuro)){
  if (grepl(',',l_platekeys_probands_notNeuro[i])){
    list_platekeys = strsplit(l_platekeys_probands_notNeuro[i], ",")[[1]]
    list_platekeys = gsub(" ", "", list_platekeys, fixed = TRUE)
    l_platekeys_probands_notNeuro_unique = c(l_platekeys_probands_notNeuro_unique,
                                             max(list_platekeys))
  }else{
    l_platekeys_probands_notNeuro_unique = c(l_platekeys_probands_notNeuro_unique,
                                             l_platekeys_probands_notNeuro[i])
  }
}
length(l_platekeys_probands_notNeuro_unique)
# 37701

# 1. Merge GRCh37 and GRCh38 info, since chromosome names are different
# GRCh38 are chr1, chr2, chr3 while GRCh37 are 1,2,3
# In SHARP everything should be GRCh38, because Andy sent us coordinates only in GRCh38
merged_data$chr = recode(merged_data$chr,
                         "1" = "chr1",
                         "2" = "chr2",
                         "3" = "chr3",
                         "4" = "chr4",
                         "5" = "chr5",
                         "6" = "chr6",
                         "7" = "chr7",
                         "8" = "chr8",
                         "9" = "chr9",
                         "10" = "chr10",
                         "11" = "chr11",
                         "12" = "chr12",
                         "13" = "chr13",
                         "14" = "chr14",
                         "15" = "chr15",
                         "16" = "chr16",
                         "17" = "chr17",
                         "18" = "chr18",
                         "19" = "chr19",
                         "20" = "chr20",
                         "21" = "chr21",
                         "22" = "chr22",
                         "X" = "chrX")

# Let's focus on SHARP genes
l_genes = unique(merged_data$gene)
length(l_genes)
# 329

# Let's focus only on SHARP genes
l_sharp = l_genes[which(grepl("SHARP", l_genes, ignore.case = TRUE))]
length(l_sharp)
# 55

sharp_merged_data = merged_data %>%
  filter(gene %in% l_sharp)
dim(sharp_merged_data)
# 5974  12

# Output folder
output_folder = 'EHv322_batch_august2020'
#dir.create(output_folder)

l_genes = unique(sharp_merged_data$gene)
for(i in 1:length(l_genes)){
  # Let's create 2 df: sharp_merged_data for ONLY PROBANDS and for ONLY PROBANDS NOT NEURO (from the list of platekeys) 
  
  sharp_merged_data_locus = sharp_merged_data %>%
    filter(gene %in% l_genes[i])
  
  sharp_merged_data_probands = data.frame()
  sharp_merged_data_probands_notNeuro = data.frame()
  for(z in 1:length(sharp_merged_data_locus$allele)){
    list_samples = strsplit(sharp_merged_data_locus$list_samples[z], ";")[[1]]
    list_samples = gsub("^EH_", "", list_samples)
    list_samples = gsub(".vcf", "", list_samples)
    list_samples = gsub(" ", "", list_samples)
    
    # Keep only in Probands and ProbandsNotNeuro
    list_samples_probands = intersect(list_samples, l_platekeys_probands_unique)
    list_samples_probands_notNeuro = intersect(list_samples, l_platekeys_probands_notNeuro_unique)
    
    new_line_probands = cbind(sharp_merged_data_locus$gene[z], sharp_merged_data_locus$allele[z], length(unique(list_samples_probands)))
    new_line_probands_notNeuro = cbind(sharp_merged_data_locus$gene[z], sharp_merged_data_locus$allele[z], length(unique(list_samples_probands_notNeuro)))
    
    sharp_merged_data_probands = rbind(sharp_merged_data_probands,
                                       new_line_probands)
    
    sharp_merged_data_probands_notNeuro = rbind(sharp_merged_data_probands_notNeuro,
                                                new_line_probands_notNeuro)
    
  }
  colnames(sharp_merged_data_probands) = c("gene", "repeat-size", "total_num_samples")
  colnames(sharp_merged_data_probands_notNeuro) = c("gene", "repeat-size", "total_num_samples")
  
  sharp_boxplot_probands = data.frame()
  for(j in 1:length(sharp_merged_data_probands$gene)){
    allele = sharp_merged_data_probands$`repeat-size`[j]
    repeat_size_probands = sharp_merged_data_probands$total_num_samples[j]
    new_line = c(l_genes[i], allele)
    sharp_boxplot_probands = rbind(sharp_boxplot_probands,
                                   as.data.frame(do.call("rbind", replicate(repeat_size_probands, new_line, simplify = FALSE))))
  }
  sharp_boxplot_probands$cohort = rep("only probands", length(sharp_boxplot_probands$V1))
  
  sharp_boxplot_probands_notNeuro = data.frame()
  for(j in 1:length(sharp_merged_data_probands_notNeuro$gene)){
    allele = sharp_merged_data_probands_notNeuro$`repeat-size`[j]
    repeat_size_probands = sharp_merged_data_probands_notNeuro$total_num_samples[j]
    new_line = c(l_genes[i], allele)
    sharp_boxplot_probands_notNeuro = rbind(sharp_boxplot_probands_notNeuro,
                                            as.data.frame(do.call("rbind", replicate(repeat_size_probands, new_line, simplify = FALSE))))
  }
  sharp_boxplot_probands_notNeuro$cohort = rep("only probands not neurology", length(sharp_boxplot_probands_notNeuro$V1))
  
  colnames(sharp_boxplot_probands) = c("gene", "repeat-size", "cohort")
  colnames(sharp_boxplot_probands_notNeuro) = c("gene", "repeat-size", "cohort")
  
  merged_sharp_boxplot = rbind(sharp_boxplot_probands,
                               sharp_boxplot_probands_notNeuro)
  
  # Create histogram for the gene
  l_gene_repeat_size = aux_df$allele
  sharp_barplot = data.frame(number_repeats = l_gene_repeat_size, af = aux_df$total_num_samples)
  
  # order by 'number of repetition'
  sharp_barplot = unique(sharp_barplot[order(sharp_barplot[,1]),])
  
  rownames(sharp_barplot) = sharp_barplot$number_repeats
  
  sharp_barplot$number_repeats = as.numeric(sharp_barplot$number_repeats)
  
  png_name = paste(l_genes[i], 'png', sep = ".")
  png_name = paste(output_folder, png_name, sep = "/")
  
  min_value = min(sharp_barplot$number_repeats)
  max_value = max(sharp_barplot$number_repeats)
  
  gene_histo = ggplot(unique(sharp_barplot), aes(x = number_repeats, y = af)) + 
    scale_x_continuous(limits=c(min_value, max_value)) +
    geom_bar(stat = "identity") + 
    geom_text(aes(label=af), vjust=0, size = 4, colour = "grey") +
    theme(legend.position = "none") +
    theme(axis.title.x=element_blank(),
          axis.title.y=element_blank(),
          axis.ticks.x=element_blank(),
          axis.ticks.y=element_blank())
    #ylab("Allele frequency") + 
    #xlab("Repeat sizes (repeat units)") + 
    #ggtitle(l_genes[i]) 
    #coord_cartesian(xlim = c(min_value,max_value))
  
  # Create boxplots for PROBANDS only and PROBANDS NOT IN NEURO for each gene
  merged_sharp_boxplot$gene = as.character(merged_sharp_boxplot$gene)
  merged_sharp_boxplot$`repeat-size` = as.integer(as.character(merged_sharp_boxplot$`repeat-size`))
  merged_sharp_boxplot$cohort = as.character(merged_sharp_boxplot$cohort)
  gene_boxplot = ggplot(merged_sharp_boxplot, aes(x = `repeat-size`, y = l_genes[i], fill = cohort)) +
    #scale_y_discrete(limits=c(min_value, max_value)) +
    geom_violin() +
    xlab("Repeat size") +
    #coord_cartesian(xlim = c(min_value,max_value)) +
    geom_boxplot(width=0.1) +
    theme(legend.position = "none") +
    theme(axis.title.x=element_blank(),
          axis.title.y=element_blank(),
          axis.ticks.x=element_blank(),
          axis.ticks.y=element_blank())
  
  # Combining histo and boxplots
  together_plot_locus = cowplot::plot_grid(gene_histo,
                                           gene_boxplot, 
                                           ncol = 1, 
                                           rel_heights = c(2, 1),
                                           align = 'v',
                                           axis = 'lr')
  png(png_name)
  print(together_plot_locus)
  dev.off()
  
}
#
merged_data = merged_data %>%
  group_by(chr, gene, allele) %>%
  mutate(total_num_samples = sum(num_samples)) %>%
  ungroup() %>%
  as.data.frame() 





merged_data_simpl = merged_data %>% 
  select(chr, gene, allele, total_num_samples)
merged_data_simpl = unique(merged_data_simpl)
dim(merged_data_simpl)
# 21013  4

