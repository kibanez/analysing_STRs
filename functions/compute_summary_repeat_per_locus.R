# Function that generates quantitative  repeat-size summary for each locus and continental group
compute_summary_repeat_per_locus <- function(df_input, gene_name, output_folder) {
  # Focus on the gene_name
  df_input = df_input %>%
    filter(gene %in% gene_name)
  
  # Create phiction dataframe in order to call a simple summary with the repeat sizes
  df_aux = data.frame()
  l_popu = unique(df_input$population)
  for (i in l_popu){
    df_number_repeats = df_input %>% 
      filter(population %in% i) %>%
      select(allele, num_samples)
    
    summary_population = rep(df_number_repeats$allele, df_number_repeats$num_samples)
    population_popu = i
    median_popu = as.integer(summary(summary_population)[3])
    q1_popu = as.integer(summary(summary_population)[2])
    q3_popu = as.integer(summary(summary_population)[5])
    
    df_aux = rbind(df_aux,
                   cbind(population_popu, median_popu, q1_popu, q3_popu))
  }
  
  # Write into a file
  output_file = paste(gene_name, "summary_quantitative_repeatsize_distribution.tsv", sep = "_")
  output_file = paste(output_folder, output_file, sep = "")
  write.table(df_aux,
              output_file,
              quote = F,
              row.names = F,
              col.names = T,
              sep = "\t")
}