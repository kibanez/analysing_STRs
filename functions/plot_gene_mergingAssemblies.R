# Function that takes already merged GRCh37 and GRCh38 tables and then plots the STR repeat-sizes distribution for a gene
plot_gene_mergingAssemblies <- function(df_input, gene_name, output_folder) {
  df_gene = df_input %>% filter(gene %in% gene_name)
  alt_number = df_gene$allele
  
  df_gene_barplot = data.frame(number_repeats = alt_number, af = df_gene$total_num_samples)
  
  # order by 'number of repetition'
  df_gene_barplot = unique(df_gene_barplot[order(df_gene_barplot[,1]),])
  
  rownames(df_gene_barplot) = df_gene_barplot$number_repeats
  
  df_gene_barplot$number_repeats = as.numeric(df_gene_barplot$number_repeats)
  
  gene_name = paste(gene_name, "merged_GRCh37_GRCh38", sep = '_')
  png_name = paste(gene_name, 'png', sep = ".")
  png_name = paste(output_folder, png_name, sep = "/")
  
  min_value = min(df_gene_barplot$number_repeats)
  max_value = max(df_gene_barplot$number_repeats)
  
  aux_plot = ggplot(unique(df_gene_barplot), aes(x = number_repeats, y = af)) + 
    geom_bar(stat = "identity") + 
    geom_text(aes(label=af), vjust=0, size = 4, colour = "grey") +
    ylab("Allele frequency") + 
    xlab("Repeat sizes (repeat units)") + 
    ggtitle(gene_name) + 
    #geom_vline(xintercept = threshold_normal, colour = 'blue', lty = 2) + 
    #geom_vline(xintercept = threshold_pathogenic, colour = 'red', lty = 2) + 
    coord_cartesian(xlim = c(min_value,max_value))
  
  png(png_name)
  print(aux_plot, res=300)
  dev.off()
}
