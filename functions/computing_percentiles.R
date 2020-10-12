# Function that given a dataframe with different genes, it computes the percentiles for each gene across given population
computing_percentiles <- function(df_input){
  l_prob = c(0,0.25,0.5,0.75,0.999,1)
  p_names <- map_chr(l_prob, ~paste0(.x*100, "%"))
  p_funs <- map(l_prob, 
                ~partial(quantile, probs = .x, na.rm = TRUE)) %>% 
    set_names(nm = p_names)
  
  df_percentiles = df_input %>%
    group_by(gene) %>%
    summarize_at(vars(allele), funs(!!!p_funs)) %>%
    ungroup() %>%
    as.data.frame()
  
 return(df_percentiles)
}