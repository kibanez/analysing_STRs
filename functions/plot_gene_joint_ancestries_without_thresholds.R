# Function that plots jointly all STR distribution across all ancestries given
# Ignoring thresholds or cutoffs
plot_gene_joint_ancestries_without_cutoff <- function(df_input, gene_name, output_folder) {
  df_gene = df_input %>% filter(gene %in% gene_name)
  alt_number = df_gene$repeat.size
  population = df_gene$population
  df_gene_barplot = data.frame(number_repeats = alt_number, af = df_gene$num_samples, population = population)
  
  pdf_name = paste(output_folder, gene_name, sep = "/")
  pdf_name = paste(pdf_name, "joint_ancestries", sep = "_")
  png_name = paste(pdf_name, 'png', sep = ".")
  pdf_name = paste(pdf_name, 'pdf', sep = ".")
  
  min_value = min(df_gene_barplot$number_repeats)
  max_value = max(df_gene_barplot$number_repeats)
  
  
  joint_plot = ggplot(df_gene_barplot, aes(x = number_repeats, y = af, group = population, color = population)) + 
    geom_line() + 
    ylab("Allele frequency") + 
    xlab("Repeat sizes (repeat units)") + 
    ggtitle(gene_name) + 
    coord_cartesian(xlim = c(min_value,max_value))
  
  pdf(pdf_name)
  print(joint_plot)	
  dev.off()
  
  png(png_name)
  print(joint_plot)
  dev.off()
  
}
