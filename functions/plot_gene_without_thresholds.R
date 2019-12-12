# Function that plots the STR repeat-size frequencies for a gene/locus across the cohort
# Not taking into account cut-offs/thresholds
plot_gene_without_cutoff <- function(df_input, gene_name, output_folder, assembly, ancestry) {
  df_gene = df_input %>% filter(gene %in% gene_name)
  
  alt_number = df_gene$allele
  
  df_gene_barplot = data.frame(number_repeats = alt_number, af = df_gene$num_samples)
  
  # order by 'number of repetition'
  df_gene_barplot = unique(df_gene_barplot[order(df_gene_barplot[,1]),])
  
  rownames(df_gene_barplot) = df_gene_barplot$number_repeats
  
  df_gene_barplot$number_repeats = as.numeric(df_gene_barplot$number_repeats)
  
  gene_name = paste(gene_name, assembly, sep = '_')
  pdf_name = paste(output_folder, gene_name, sep = "/")
  pdf_name = paste(pdf_name, ancestry, sep = "_")
  png_name = paste(pdf_name, 'png', sep = ".")
  pdf_name = paste(pdf_name, 'pdf', sep = ".")
  
  min_value = min(df_gene_barplot$number_repeats)
  max_value = max(threshold_pathogenic + 1, df_gene_barplot$number_repeats)
  
  aux_plot = ggplot(unique(df_gene_barplot), aes(x = number_repeats, y = af)) + 
    geom_bar(stat = "identity") + 
    geom_text(aes(label=af), vjust=0, size = 4, colour = "grey") +
    ylab("Allele frequency") + 
    xlab("Repeat sizes (repeat units)") + 
    ggtitle(gene_name) + 
    coord_cartesian(xlim = c(min_value,max_value))
  
  pdf(pdf_name)
  print(aux_plot)	
  dev.off()
  
  png(png_name)
  print(aux_plot)
  dev.off()
}
