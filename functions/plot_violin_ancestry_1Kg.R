plot_violin_ancestry_1Kg <- function(df_input, gene_name, gene_data_normal, gene_data_pathogenic, output_folder){
  threshold_normal = gene_data_normal %>% filter(grepl(gene_name, locus)) %>% select(threshold) %>% unlist() %>% unname()
  threshold_pathogenic = gene_data_pathogenic %>% filter(grepl(gene_name, locus)) %>% select(threshold) %>% unlist() %>% unname()

  df_gene = df_input %>% 
    filter(gene %in% gene_name) %>%
    select(allele, num_samples, population, superpopulation)
  
  df_gene2 = data.frame()
  l_ancestries = unique(df_gene$superpopulation)
  for(i in 1:length(l_ancestries)){
    aux_df = df_gene %>% filter(superpopulation %in% l_ancestries[i])
    for(j in 1:length(aux_df$allele)){
      allele = aux_df$allele[j]
      repeat_size = aux_df$num_samples[j]
      superpopulation = aux_df$superpopulation[j]
      subpopulation = aux_df$population[j]
      
      new_line = c(allele, subpopulation, superpopulation)
      
      df_gene2 = rbind(df_gene2,
                       as.data.frame(do.call("rbind", replicate(repeat_size, new_line, simplify = FALSE))))
    }
  }
  
  colnames(df_gene2) = c("repeat_size", "subpopulation" , "superpopulation")
  df_gene2$repeat_size = as.integer(as.character(df_gene2$repeat_size))
  df_gene2$subpopulation = as.character(df_gene2$subpopulation)
  df_gene2$superpopulation = as.character(df_gene2$superpopulation)
  
  pdf_name = paste(output_folder, gene_name, sep = "/")
  pdf_name = paste(pdf_name, "joint_ancestries_violin_plot", sep = "_")
  png_name = paste(pdf_name, "png", sep = ".")
  pdf_name = paste(pdf_name, 'pdf', sep = ".")
  
  violin_plot = ggplot(df_gene2, aes(x = subpopulation, y=repeat_size, fill = superpopulation)) +
    geom_violin() +
    coord_flip() +
    xlab("Sub-population") + 
    ylab("Repeat sizes (repeat units)") + 
    ggtitle(gene_name) +
    geom_boxplot(width=0.1) +
    guides(fill=guide_legend(reverse=TRUE)) +
    scale_x_discrete(limits=c("ACB", "ASW","ESN", "GWD", "LWK", "MSL", "YRI",
                              "CLM", "MXL", "PEL","PUR",
                              "CDX", "CHB", "CHS", "JPT", "KHV",
                              "CEU", "GBR", "FIN", "IBS", "TSI", 
                              "BEB", "GIH", "ITU", "PJL", "STU"))
  
  pdf(pdf_name)
  print(violin_plot)
  dev.off()
  
  png(png_name)
  print(violin_plot)
  dev.off()
}
