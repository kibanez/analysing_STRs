# Function that plots jointly all STR distribution across all ancestries given
plot_gene_joint_ancestries_1Kg <- function(df_input, gene_name, gene_data_normal, gene_data_pathogenic, output_folder) {
  threshold_normal = gene_data_normal %>% filter(grepl(gene_name, locus)) %>% select(threshold) %>% unlist() %>% unname()
  threshold_pathogenic = gene_data_pathogenic %>% filter(grepl(gene_name, locus)) %>% select(threshold) %>% unlist() %>% unname()
  
  df_gene = df_input %>% filter(gene %in% gene_name)
  alt_number = df_gene$allele
  subpopulation = df_gene$population
  superpopulation = df_gene$superpopulation
  df_gene_barplot = data.frame(number_repeats = alt_number, af = df_gene$num_samples, 
                               subpopulation = subpopulation, superpopulation = superpopulation)
  
  pdf_name = paste(output_folder, gene_name, sep = "/")
  pdf_name = paste(pdf_name, "joint_ancestries", sep = "_")
  png_name = paste(pdf_name, 'png', sep = ".")
  pdf_name = paste(pdf_name, 'pdf', sep = ".")
  
  min_value = min(df_gene_barplot$number_repeats)
  max_value = max(threshold_pathogenic + 1, df_gene_barplot$number_repeats)
  
  
  joint_plot = ggplot(df_gene_barplot, aes(x = number_repeats, y = af, group = subpopulation, color = superpopulation)) + 
    geom_line() + 
    ylab("Allele frequency") + 
    xlab("Repeat sizes (repeat units)") + 
    ggtitle(gene_name) + 
    #geom_vline(xintercept = threshold_normal, colour = 'blue', lty = 2) + 
    #geom_vline(xintercept = threshold_pathogenic, colour = 'red', lty = 2) + 
    coord_cartesian(xlim = c(min_value,max_value))
  
  pdf(pdf_name)
  print(joint_plot)	
  dev.off()
  
  png(png_name)
  print(joint_plot)
  dev.off()
  
}
