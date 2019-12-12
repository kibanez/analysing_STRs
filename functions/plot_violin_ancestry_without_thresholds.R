# Adapt plot_violin_ancestry function not taking cut-offs
plot_violin_ancestry_without_cutoff <- function(df_input, gene_name, output_folder){
  df_gene = df_input %>% 
    filter(gene %in% gene_name) %>%
    select(allele, num_samples, population)
  
  df_gene2 = data.frame()
  l_ancestries = unique(df_gene$population)
  for(i in 1:length(l_ancestries)){
    aux_df = df_gene %>% filter(population %in% l_ancestries[i])
    for(j in 1:length(aux_df$allele)){
      allele = aux_df$allele[j]
      repeat_size = aux_df$num_samples[j]
      population = aux_df$population[j]
      
      new_line = c(allele, population)
      
      df_gene2 = rbind(df_gene2,
                       as.data.frame(do.call("rbind", replicate(repeat_size, new_line, simplify = FALSE))))
    }
  }
  
  dim(df_gene2)
  
  colnames(df_gene2) = c("repeat_size", "population")
  df_gene2$repeat_size = as.integer(as.character(df_gene2$repeat_size))
  df_gene2$population = as.character(df_gene2$population)
  
  pdf_name = paste(output_folder, gene_name, sep = "/")
  pdf_name = paste(pdf_name, "joint_ancestries_violin_plot", sep = "_")
  png_name = paste(pdf_name, "png", sep = ".")
  pdf_name = paste(pdf_name, 'pdf', sep = ".")
  
  my_comparisons=list(c("AFR","AMR"), c("AFR","EUR"), c("AFR","EAS"), c("AFR","ASI"))
  #                    c("AMR","EUR"), c("AMR","EAS"), c("AMR","ASI"),
  #                    c("EUR", "EAS"), c("EUR", "ASI"),
  #                    c("EAS", "ASI"))
  violin_plot = ggplot(df_gene2, aes(x = population, y=repeat_size, fill = population)) +
    geom_violin() +
    coord_flip() +
    xlab("Population") + 
    ylab("Repeat sizes (repeat units)") + 
    ggtitle(gene_name) 
    #stat_compare_means(comparisons = my_comparisons,
#                       method = "wilcox.test")
  
  pdf(pdf_name)
  print(violin_plot)
  dev.off()
  
  png(png_name)
  print(violin_plot)
  dev.off()
}
