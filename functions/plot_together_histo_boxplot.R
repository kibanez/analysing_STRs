# Function that plots histogram together with boxplot
plot_together_histo_boxplot <- function(df_input, gene_name, output_folder, l_platekeys_probands_neuro_unique, l_platekeys_probands_neuro_notNeuro_unique) {
  # Let's create 2 df: sharp_merged_data for ONLY PROBANDS and for ONLY PROBANDS NOT NEURO (from the list of platekeys) 
  df_locus = df_input %>%
    filter(gene %in% gene_name)
  
  df_all_genomes = df_locus
  # Since chr are being merged, we need to compute total_num_samples
  df_all_genomes = df_all_genomes %>%
    group_by(allele) %>%
    mutate(total_num_samples = sum(num_samples)) %>%
    ungroup() %>%
    select(allele, total_num_samples) %>%
    as.data.frame() %>%
    unique()
  
  # merge chr and 
  df_locus = df_locus %>%
    group_by(chr, gene, allele) %>%
    summarise(list_samples_merged = paste(list_samples, collapse=";")) %>%
    ungroup() %>%
    as.data.frame() 
  
  df_locus = df_locus %>%
    select(gene, allele, list_samples_merged)
  df_locus = unique(df_locus)
  
  df_probands = data.frame()
  df_probands_notNeuro = data.frame()
  for(z in 1:length(df_locus$allele)){
    list_samples = strsplit(df_locus$list_samples_merged[z], ";")[[1]]
    index_x2 = which(grepl("x2", list_samples))
    
    list_samples = gsub("^EH_", "", list_samples)
    list_samples = gsub(".vcf", "", list_samples)
    list_samples = gsub("_x2", "", list_samples)
    list_samples = gsub(" ", "", list_samples)
    
    # Keep only in Probands and ProbandsNotNeuro
    list_samples_probands = unique(intersect(list_samples, l_platekeys_probands_neuro_unique))
    list_samples_probands_notNeuro = unique(intersect(list_samples, l_platekeys_probands_notNeuro_unique))
    
    list_samples_probands_x2 = intersect(list_samples[index_x2], l_platekeys_probands_neuro_unique)
    list_samples_probands_notNeuro_x2 = intersect(list_samples[index_x2], l_platekeys_probands_notNeuro_unique)
    
    new_line_probands = cbind(df_locus$gene[z], df_locus$allele[z], length(c(list_samples_probands, list_samples_probands_x2)))
    new_line_probands_notNeuro = cbind(df_locus$gene[z], df_locus$allele[z], length(c(list_samples_probands_notNeuro, list_samples_probands_notNeuro_x2)))
    
    df_probands = rbind(df_probands,
                        new_line_probands)
    
    df_probands_notNeuro = rbind(df_probands_notNeuro,
                                 new_line_probands_notNeuro)
    
  }
  colnames(df_probands) = c("gene", "repeat_size", "total_num_samples")
  colnames(df_probands_notNeuro) = c("gene", "repeat_size", "total_num_samples")
  
  # Filter out those repeat-sizes for which `total_num_samples` == 0 (after filtering out by probands and probands not neuro)
  df_probands = df_probands %>%
    filter(total_num_samples != 0)

  df_probands_notNeuro = df_probands_notNeuro %>%
    filter(total_num_samples != 0)
  
  df_boxplot_probands = data.frame()
  for(j in 1:length(df_probands$gene)){
    allele = as.integer(as.character(df_probands$repeat_size[j]))
    repeat_size_probands = as.integer(as.character(df_probands$total_num_samples[j]))
    new_line = c(l_genes[i], allele)
    df_boxplot_probands = rbind(df_boxplot_probands,
                                   as.data.frame(do.call("rbind", replicate(repeat_size_probands, new_line, simplify = FALSE))))
  }
  df_boxplot_probands$cohort = rep("only probands in neurology", length(df_boxplot_probands$V1))
  
  df_boxplot_probands_notNeuro = data.frame()
  for(j in 1:length(df_probands_notNeuro$gene)){
    allele = as.integer(as.character(df_probands_notNeuro$repeat_size[j]))
    repeat_size_probands = as.integer(as.character(df_probands_notNeuro$total_num_samples[j]))
    new_line = c(l_genes[i], allele)
    df_boxplot_probands_notNeuro = rbind(df_boxplot_probands_notNeuro,
                                            as.data.frame(do.call("rbind", replicate(repeat_size_probands, new_line, simplify = FALSE))))
  }
  df_boxplot_probands_notNeuro$cohort = rep("only probands not neurology", length(df_boxplot_probands_notNeuro$V1))
  
  colnames(df_boxplot_probands) = c("gene", "repeat_size", "cohort")
  colnames(df_boxplot_probands_notNeuro) = c("gene", "repeat_size", "cohort")
  
  merged_df_boxplot = rbind(df_boxplot_probands,
                               df_boxplot_probands_notNeuro)
  
  # Create boxplots for PROBANDS only and PROBANDS NOT IN NEURO for each gene
  merged_df_boxplot$gene = as.character(merged_df_boxplot$gene)
  merged_df_boxplot$repeat_size = as.integer(as.character(merged_df_boxplot$repeat_size))
  merged_df_boxplot$cohort = as.character(merged_df_boxplot$cohort)
  
  gene_boxplot = ggplot(merged_df_boxplot, aes(x = repeat_size, y = cohort, fill = cohort)) +
    geom_violin() +
    xlab("Repeat size") +
    geom_boxplot(width=0.1) +
    theme(legend.position = "none") +
    theme(axis.title.x=element_blank(),
          axis.title.y=element_blank(),
          axis.ticks.x=element_blank(),
          axis.ticks.y=element_blank())
  
  gene_boxplot2 = merged_df_boxplot %>% 
    ggplot(aes(x=repeat_size,y=cohort, fill=cohort)) +
    geom_boxplot() + geom_jitter(width=0.1,alpha=0.2) +
    theme(legend.position = "none") +
    theme(axis.title.x=element_blank(),
          axis.title.y=element_blank(),
          axis.ticks.x=element_blank(),
          axis.ticks.y=element_blank())
  
  # Create histogram for the gene
  l_gene_repeat_size_probands = as.integer(as.character(df_probands$repeat_size))
  l_gene_repeat_size_probands_notNeuro = as.integer(as.character(df_probands_notNeuro$repeat_size))
  genes_barplot_probands = data.frame(number_repeats = l_gene_repeat_size_probands,
                                      af = df_probands$total_num_samples)
  
  genes_barplot_probands_notNeuro = data.frame(number_repeats = l_gene_repeat_size_probands_notNeuro,
                                               af = df_probands_notNeuro$total_num_samples)
  
  # order by 'number of repetition'
  genes_barplot_probands = unique(genes_barplot_probands[order(genes_barplot_probands[,1]),])
  genes_barplot_probands_notNeuro = unique(genes_barplot_probands_notNeuro[order(genes_barplot_probands_notNeuro[,1]),])
  
  rownames(genes_barplot_probands) = genes_barplot_probands$number_repeats
  rownames(genes_barplot_probands_notNeuro) = genes_barplot_probands_notNeuro$number_repeats
  
  genes_barplot_probands$number_repeats = as.numeric(genes_barplot_probands$number_repeats)
  genes_barplot_probands_notNeuro$number_repeats = as.numeric(genes_barplot_probands_notNeuro$number_repeats)
  
  genes_barplot_probands$cohort = rep("only probands", length(genes_barplot_probands$number_repeats))
  genes_barplot_probands_notNeuro$cohort = rep("only probands not neurology", length(genes_barplot_probands_notNeuro$number_repeats))
  
  sharp_barplot_merged = rbind(genes_barplot_probands,
                               genes_barplot_probands_notNeuro)
  
  sharp_barplot_merged$af = as.numeric(as.character(sharp_barplot_merged$af))
  
  png_name = paste(gene_name, 'png', sep = ".")
  png_name = paste(output_folder, png_name, sep = "/")
  
  min_value = min(genes_barplot_probands$number_repeats, genes_barplot_probands_notNeuro$number_repeats)
  max_value = max(genes_barplot_probands$number_repeats, genes_barplot_probands_notNeuro$number_repeats)
  
  # Histogram with 93k genomes (all genomes, not distinguising them between probands neuro, probands not in neuro, etc.)
  df_all_genomes_barplot = data.frame(number_repeats = df_all_genomes$allele, 
                                      af = df_all_genomes$total_num_samples)
  
  # order by 'number of repetition'
  df_all_genomes_barplot = unique(df_all_genomes[order(df_all_genomes_barplot[,1]),])
  rownames(df_all_genomes_barplot) = df_all_genomes_barplot$number_repeats
  min_value = min(df_all_genomes_barplot$allele)
  max_value = max(df_all_genomes_barplot$allele)
  
  gene_histo = ggplot(df_all_genomes_barplot, aes(x = allele, y = total_num_samples)) + 
    scale_x_continuous(limits=c(min_value, max_value)) +
    geom_bar(stat = "identity") + 
    geom_text(aes(label=total_num_samples), vjust=0, size = 4, colour = "grey") +
    theme(legend.position = "none") +
    theme(axis.title.x=element_blank(),
          axis.title.y=element_blank(),
          axis.ticks.x=element_blank(),
          axis.ticks.y=element_blank())
  
    #gene_histo = ggplot(unique(sharp_barplot_merged), aes(x = number_repeats, y = af, fill = cohort)) + 
  #  scale_x_continuous(limits=c(min_value, max_value)) +
  #  geom_bar(stat = "identity") + 
  #  geom_text(aes(label=af), vjust=0, size = 4, colour = "grey") +
  #  theme(legend.position = "none") +
  #  theme(axis.title.x=element_blank(),
  #        axis.title.y=element_blank(),
  #        axis.ticks.x=element_blank(),
  #        axis.ticks.y=element_blank())
  
  # Combining histo and boxplots
  together_plot_locus = cowplot::plot_grid(gene_histo,
                                           gene_boxplot, 
                                           ncol = 1, 
                                           rel_heights = c(2, 1),
                                           align = 'v',
                                           axis = 'lr')
  png(png_name, units="in", width=5, height=5, res=300)
  print(together_plot_locus)
  dev.off()
}
