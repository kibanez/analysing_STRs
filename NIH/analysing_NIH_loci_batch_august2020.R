# Objective: analyse loci analysed by NIH within 100kGP, after EHdn 
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
setwd("~/Documents/STRs/ANALYSIS/NIH/")

# load merged august data
merged_data = read.csv("~/Documents/STRs/data/research/batch_august2020/output_EHv3.2.2_vcfs/merged/merged_93446_genomes_EHv322_batch_august2020.tsv",
                       stringsAsFactors = F,
                       sep = "\t",
                       header = T)
dim(merged_data)
# 27238  12

# 1. Merge GRCh37 and GRCh38 info, since chromosome names are different
# GRCh38 are chr1, chr2, chr3 while GRCh37 are 1,2,3
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

# Let's focus on NIH loci - the ones starting by `^chr`
l_genes = unique(merged_data$gene)
length(l_genes)
# 329

l_nih = l_genes[which(grepl("^chr", l_genes, ignore.case = TRUE))]
length(l_nih)
# 192

# Output folder
output_folder = 'EHv322_batch_august2020'
dir.create(output_folder)

for (i in 1:length(l_sharp)){
  plot_gene_mergingAssemblies(merged_data_simpl, l_sharp[i], output_folder)
}

# Plot boxplots across all loci
sharp_merged_data = merged_data_simpl %>%
  filter(gene %in% l_sharp)
dim(sharp_merged_data)
# 3367  4

l_genes = unique(sharp_merged_data$gene)
l_genes = unique(sharp_merged_data$gene)
for(i in 1:length(l_genes)){
  # Let's create 2 df: sharp_merged_data for ONLY PROBANDS and for ONLY PROBANDS NOT NEURO (from the list of platekeys) 
  sharp_merged_data_locus = sharp_merged_data %>%
    filter(gene %in% l_genes[i])
  
  # merge chr and 
  sharp_merged_data_locus = sharp_merged_data_locus %>%
    group_by(chr, gene, allele) %>%
    summarise(list_samples_merged = paste(list_samples, collapse=";")) %>%
    ungroup() %>%
    as.data.frame() 
  
  sharp_merged_data_locus = sharp_merged_data_locus %>%
    select(gene, allele, list_samples_merged)
  sharp_merged_data_locus = unique(sharp_merged_data_locus)
  
  sharp_merged_data_probands = data.frame()
  sharp_merged_data_probands_notNeuro = data.frame()
  for(z in 1:length(sharp_merged_data_locus$allele)){
    list_samples = strsplit(sharp_merged_data_locus$list_samples_merged[z], ";")[[1]]
    index_x2 = which(grepl("x2", list_samples))
    
    list_samples = gsub("^EH_", "", list_samples)
    list_samples = gsub(".vcf", "", list_samples)
    list_samples = gsub("_x2", "", list_samples)
    list_samples = gsub(" ", "", list_samples)
    
    # Keep only in Probands and ProbandsNotNeuro
    list_samples_probands = unique(intersect(list_samples, l_platekeys_probands_unique))
    list_samples_probands_notNeuro = unique(intersect(list_samples, l_platekeys_probands_notNeuro_unique))
    
    list_samples_probands_x2 = intersect(list_samples[index_x2], l_platekeys_probands_unique)
    list_samples_probands_notNeuro_x2 = intersect(list_samples[index_x2], l_platekeys_probands_notNeuro_unique)
    
    new_line_probands = cbind(sharp_merged_data_locus$gene[z], sharp_merged_data_locus$allele[z], length(c(list_samples_probands, list_samples_probands_x2)))
    new_line_probands_notNeuro = cbind(sharp_merged_data_locus$gene[z], sharp_merged_data_locus$allele[z], length(c(list_samples_probands_notNeuro, list_samples_probands_notNeuro_x2)))
    
    sharp_merged_data_probands = rbind(sharp_merged_data_probands,
                                       new_line_probands)
    
    sharp_merged_data_probands_notNeuro = rbind(sharp_merged_data_probands_notNeuro,
                                                new_line_probands_notNeuro)
    
  }
  colnames(sharp_merged_data_probands) = c("gene", "repeat_size", "total_num_samples")
  colnames(sharp_merged_data_probands_notNeuro) = c("gene", "repeat_size", "total_num_samples")
  
  # Filter out those repeat-sizes for which `total_num_samples` == 0 (after filtering out by probands and probands not neuro)
  #sharp_merged_data_probands = sharp_merged_data_probands %>%
  #  filter(total_num_samples != 0)
  
  #sharp_merged_data_probands_notNeuro = sharp_merged_data_probands_notNeuro %>%
  #  filter(total_num_samples !=0)
  
  sharp_boxplot_probands = data.frame()
  for(j in 1:length(sharp_merged_data_probands$gene)){
    allele = as.integer(as.character(sharp_merged_data_probands$repeat_size[j]))
    repeat_size_probands = as.integer(as.character(sharp_merged_data_probands$total_num_samples[j]))
    new_line = c(l_genes[i], allele)
    sharp_boxplot_probands = rbind(sharp_boxplot_probands,
                                   as.data.frame(do.call("rbind", replicate(repeat_size_probands, new_line, simplify = FALSE))))
  }
  sharp_boxplot_probands$cohort = rep("only probands", length(sharp_boxplot_probands$V1))
  
  sharp_boxplot_probands_notNeuro = data.frame()
  for(j in 1:length(sharp_merged_data_probands_notNeuro$gene)){
    allele = as.integer(as.character(sharp_merged_data_probands_notNeuro$repeat_size[j]))
    repeat_size_probands = as.integer(as.character(sharp_merged_data_probands_notNeuro$total_num_samples[j]))
    new_line = c(l_genes[i], allele)
    sharp_boxplot_probands_notNeuro = rbind(sharp_boxplot_probands_notNeuro,
                                            as.data.frame(do.call("rbind", replicate(repeat_size_probands, new_line, simplify = FALSE))))
  }
  sharp_boxplot_probands_notNeuro$cohort = rep("only probands not neurology", length(sharp_boxplot_probands_notNeuro$V1))
  
  colnames(sharp_boxplot_probands) = c("gene", "repeat_size", "cohort")
  colnames(sharp_boxplot_probands_notNeuro) = c("gene", "repeat_size", "cohort")
  
  merged_sharp_boxplot = rbind(sharp_boxplot_probands,
                               sharp_boxplot_probands_notNeuro)
  
  # Create boxplots for PROBANDS only and PROBANDS NOT IN NEURO for each gene
  merged_sharp_boxplot$gene = as.character(merged_sharp_boxplot$gene)
  merged_sharp_boxplot$repeat_size = as.integer(as.character(merged_sharp_boxplot$repeat_size))
  merged_sharp_boxplot$cohort = as.character(merged_sharp_boxplot$cohort)
  
  gene_boxplot = ggplot(merged_sharp_boxplot, aes(x = repeat_size, y = cohort, fill = cohort)) +
    geom_violin() +
    xlab("Repeat size") +
    geom_boxplot(width=0.1) +
    theme(legend.position = "none") +
    theme(axis.title.x=element_blank(),
          axis.title.y=element_blank(),
          axis.ticks.x=element_blank(),
          axis.ticks.y=element_blank())
  
  gene_boxplot2 = merged_sharp_boxplot %>% 
    ggplot(aes(x=repeat_size,y=cohort, fill=cohort)) +
    geom_boxplot() + geom_jitter(width=0.1,alpha=0.2) +
    theme(legend.position = "none") +
    theme(axis.title.x=element_blank(),
          axis.title.y=element_blank(),
          axis.ticks.x=element_blank(),
          axis.ticks.y=element_blank())
  
  # Create histogram for the gene
  l_gene_repeat_size_probands = sharp_merged_data_probands$repeat_size
  l_gene_repeat_size_probands_notNeuro = sharp_merged_data_probands_notNeuro$repeat_size
  sharp_barplot_probands = data.frame(number_repeats = l_gene_repeat_size_probands,
                                      af = sharp_merged_data_probands$total_num_samples)
  
  sharp_barplot_probands_notNeuro = data.frame(number_repeats = l_gene_repeat_size_probands_notNeuro,
                                               af = sharp_merged_data_probands_notNeuro$total_num_samples)
  
  # order by 'number of repetition'
  sharp_barplot_probands = unique(sharp_barplot_probands[order(sharp_barplot_probands[,1]),])
  sharp_barplot_probands_notNeuro = unique(sharp_barplot_probands_notNeuro[order(sharp_barplot_probands_notNeuro[,1]),])
  
  rownames(sharp_barplot_probands) = sharp_barplot_probands$number_repeats
  rownames(sharp_barplot_probands_notNeuro) = sharp_barplot_probands_notNeuro$number_repeats
  
  sharp_barplot_probands$number_repeats = as.numeric(sharp_barplot_probands$number_repeats)
  sharp_barplot_probands_notNeuro$number_repeats = as.numeric(sharp_barplot_probands_notNeuro$number_repeats)
  
  sharp_barplot_probands$cohort = rep("only probands", length(sharp_barplot_probands$number_repeats))
  sharp_barplot_probands_notNeuro$cohort = rep("only probands not neurology", length(sharp_barplot_probands_notNeuro$number_repeats))
  
  sharp_barplot_merged = rbind(sharp_barplot_probands,
                               sharp_barplot_probands_notNeuro)
  
  sharp_barplot_merged$af = as.numeric(as.character(sharp_barplot_merged$af))
  
  png_name = paste(l_genes[i], 'png', sep = ".")
  png_name = paste(output_folder, png_name, sep = "/")
  
  min_value = min(sharp_barplot_probands$number_repeats, sharp_barplot_probands_notNeuro$number_repeats)
  max_value = max(sharp_barplot_probands$number_repeats, sharp_barplot_probands_notNeuro$number_repeats)
  
  gene_histo = ggplot(unique(sharp_barplot_merged), aes(x = number_repeats, y = af, fill = cohort)) + 
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