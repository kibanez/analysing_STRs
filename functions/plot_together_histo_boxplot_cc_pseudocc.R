# Function that plots histogram together with boxplot
plot_together_histo_boxplot_cc_pseudocc <- function(df_input, gene_name, output_folder, l_cases, l_controls, l_pseudocases, l_pseudocontrols) {
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
  
  df_cases = data.frame()
  df_controls = data.frame()
  df_pseudocases = data.frame()
  df_pseudocontrols = data.frame()
  for(z in 1:length(df_locus$allele)){
    list_samples = strsplit(df_locus$list_samples_merged[z], ";")[[1]]
    index_x2 = which(grepl("x2", list_samples))
    
    list_samples = gsub("^EH_", "", list_samples)
    list_samples = gsub(".vcf", "", list_samples)
    list_samples = gsub("_x2", "", list_samples)
    list_samples = gsub(" ", "", list_samples)
    
    # Keep only in cases/controls/pseudocases/pseudocontrols
    list_cases = unique(intersect(list_samples, l_cases))
    list_controls = unique(intersect(list_samples, l_controls))
    list_pseudocases = unique(intersect(list_samples, l_pseudocases))
    list_pseudocontrols = unique(intersect(list_samples, l_pseudocontrols))
    
    list_cases_x2 = intersect(list_samples[index_x2], l_cases)
    list_controls_x2 = intersect(list_samples[index_x2], l_controls)
    list_pseudocases_x2 = intersect(list_samples[index_x2], l_pseudocases)
    list_pseudocontrols_x2 = intersect(list_samples[index_x2], l_pseudocontrols)
    
    new_line_cases = cbind(df_locus$gene[z], df_locus$allele[z], length(c(l_cases, list_cases_x2)))
    new_line_controls = cbind(df_locus$gene[z], df_locus$allele[z], length(c(l_controls, list_controls_x2)))
    new_line_pseudocases = cbind(df_locus$gene[z], df_locus$allele[z], length(c(l_pseudocases, list_pseudocases_x2)))
    new_line_pseudocontrols = cbind(df_locus$gene[z], df_locus$allele[z], length(c(l_pseudocontrols, list_pseudocontrols_x2)))

    df_cases = rbind(df_cases,
                     new_line_cases)
    df_controls = rbind(df_controls,
                        new_line_controls)    
    df_pseudocases = rbind(df_pseudocases,
                           new_line_pseudocases)
    df_pseudocontrols = rbind(df_pseudocontrols,
                              new_line_pseudocontrols)
  }
  colnames(df_cases) = c("gene", "repeat_size", "total_num_samples")
  colnames(df_controls) = c("gene", "repeat_size", "total_num_samples")
  colnames(df_pseudocases) = c("gene", "repeat_size", "total_num_samples")
  colnames(df_pseudocontrols) = c("gene", "repeat_size", "total_num_samples")
  
  # Filter out those repeat-sizes for which `total_num_samples` == 0 (after filtering out by probands and probands not neuro)
  df_cases = df_cases %>%
    filter(total_num_samples != 0)
  df_controls = df_controls %>%
    filter(total_num_samples != 0)
  df_pseudocases = df_pseudocases %>%
    filter(total_num_samples != 0)
  df_pseudocontrols = df_pseudocontrols %>%
    filter(total_num_samples != 0)

  df_boxplot_cases = data.frame()
  for(j in 1:length(df_cases$gene)){
    allele = as.integer(as.character(df_cases$repeat_size[j]))
    repeat_size_cases = as.integer(as.character(df_cases$total_num_samples[j]))
    new_line = c(l_genes[i], allele)
    df_boxplot_cases = rbind(df_boxplot_cases,
                             as.data.frame(do.call("rbind", replicate(repeat_size_cases, new_line, simplify = FALSE))))
  }
  df_boxplot_cases$cohort = rep("cases", length(df_boxplot_cases$V1))
  
  df_boxplot_controls = data.frame()
  for(j in 1:length(df_controls$gene)){
    allele = as.integer(as.character(df_controls$repeat_size[j]))
    repeat_size_controls = as.integer(as.character(df_controls$total_num_samples[j]))
    new_line = c(l_genes[i], allele)
    df_boxplot_controls = rbind(df_boxplot_controls,
                                as.data.frame(do.call("rbind", replicate(repeat_size_controls, new_line, simplify = FALSE))))
  }
  df_boxplot_controls$cohort = rep("controls", length(df_boxplot_controls$V1))

  df_boxplot_pseudocases = data.frame()
  for(j in 1:length(df_pseudocases$gene)){
    allele = as.integer(as.character(df_pseudocases$repeat_size[j]))
    repeat_size_pseudocases = as.integer(as.character(df_pseudocases$total_num_samples[j]))
    new_line = c(l_genes[i], allele)
    df_boxplot_pseudocases = rbind(df_boxplot_pseudocases,
                             as.data.frame(do.call("rbind", replicate(repeat_size_pseudocases, new_line, simplify = FALSE))))
  }
  df_boxplot_pseudocases$cohort = rep("pseudocases", length(df_boxplot_pseudocases$V1))

  df_boxplot_pseudocontrols = data.frame()
  for(j in 1:length(df_pseudocontrols$gene)){
    allele = as.integer(as.character(df_pseudocontrols$repeat_size[j]))
    repeat_size_pseudocontrols = as.integer(as.character(df_pseudocontrols$total_num_samples[j]))
    new_line = c(l_genes[i], allele)
    df_boxplot_pseudocontrols = rbind(df_boxplot_pseudocontrols,
                                      as.data.frame(do.call("rbind", replicate(repeat_size_pseudocontrols, new_line, simplify = FALSE))))
  }
  df_boxplot_pseudocontrols$cohort = rep("pseudocontrol", length(df_boxplot_pseudocontrols$V1))
  
  colnames(df_boxplot_cases) = c("gene", "repeat_size", "cohort")
  colnames(df_boxplot_controls) = c("gene", "repeat_size", "cohort")
  colnames(df_boxplot_pseudocases) = c("gene", "repeat_size", "cohort")
  colnames(df_boxplot_pseudocontrols) = c("gene", "repeat_size", "cohort")
  
  merged_df_boxplot = rbind(df_boxplot_cases,
                            df_boxplot_controls,
                            df_boxplot_pseudocases,
                            df_boxplot_pseudocontrols)
  
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
  l_gene_repeat_size_cases = as.integer(as.character(df_cases$repeat_size))
  l_gene_repeat_size_controls = as.integer(as.character(df_controls$repeat_size))
  l_gene_repeat_size_pseudocases = as.integer(as.character(df_pseudocases$repeat_size))
  l_gene_repeat_size_pseudocontrols = as.integer(as.character(df_pseudocontrols$repeat_size))
  
  genes_barplot_cases = data.frame(number_repeats = l_gene_repeat_size_cases,
                                   af = df_cases$total_num_samples)
  genes_barplot_controls = data.frame(number_repeats = l_gene_repeat_size_controls,
                                      af = df_controls$total_num_samples)
  genes_barplot_pseudocases = data.frame(number_repeats = l_gene_repeat_size_pseudocases,
                                         af = df_pseudocases$total_num_samples)
  genes_barplot_pseudocontrols = data.frame(number_repeats = l_gene_repeat_size_pseudocontrols,
                                            af = df_pseudocontrols$total_num_samples)
  
  # order by 'number of repetition'
  genes_barplot_cases = unique(genes_barplot_cases[order(genes_barplot_cases[,1]),])
  genes_barplot_controls = unique(genes_barplot_controls[order(genes_barplot_controls[,1]),])
  genes_barplot_pseudocases = unique(genes_barplot_pseudocases[order(genes_barplot_pseudocases[,1]),])
  genes_barplot_pseudocontrols = unique(genes_barplot_pseudocontrols[order(genes_barplot_pseudocontrols[,1]),])
  
  rownames(genes_barplot_cases) = genes_barplot_cases$number_repeats
  rownames(genes_barplot_controls) = genes_barplot_controls$number_repeats
  rownames(genes_barplot_pseudocases) = genes_barplot_pseudocases$number_repeats
  rownames(genes_barplot_pseudocontrols) = genes_barplot_pseudocontrols$number_repeats
  
  genes_barplot_cases$number_repeats = as.numeric(genes_barplot_cases$number_repeats)
  genes_barplot_controls$number_repeats = as.numeric(genes_barplot_controls$number_repeats)
  genes_barplot_pseudocases$number_repeats = as.numeric(genes_barplot_pseudocases$number_repeats)
  genes_barplot_pseudocontrols$number_repeats = as.numeric(genes_barplot_pseudocontrols$number_repeats)
  
  genes_barplot_cases$cohort = rep("cases", length(genes_barplot_cases$number_repeats))
  genes_barplot_controls$cohort = rep("controls", length(genes_barplot_controls$number_repeats))
  genes_barplot_pseudocases$cohort = rep("pseudocases", length(genes_barplot_pseudocases$number_repeats))
  genes_barplot_pseudocontrols$cohort = rep("pseudocontrols", length(genes_barplot_pseudocontrols$number_repeats))
  
  sharp_barplot_merged = rbind(genes_barplot_cases,
                               genes_barplot_controls,
                               genes_barplot_pseudocases,
                               genes_barplot_pseudocontrols)
  
  sharp_barplot_merged$af = as.numeric(as.character(sharp_barplot_merged$af))
  
  png_name = paste(gene_name, 'png', sep = ".")
  png_name = paste(output_folder, png_name, sep = "/")
  
  min_value = min(genes_barplot_cases$number_repeats, genes_barplot_controls$number_repeats, genes_barplot_pseudocases$number_repeats, genes_barplot_pseudocontrols$number_repeats)
  max_value = max(genes_barplot_cases$number_repeats, genes_barplot_controls$number_repeats, genes_barplot_pseudocases$number_repeats, genes_barplot_pseudocontrols$number_repeats)
  
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
    #geom_text(aes(label=total_num_samples), vjust=0, size = 4, colour = "grey") +
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
  png(png_name, units="in", width=5, height=5, res=300)
  print(together_plot_locus)
  dev.off()
}
