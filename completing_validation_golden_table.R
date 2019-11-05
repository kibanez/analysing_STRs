# Objective: Complete our validation table together with the max CI value estimated by EH for each allele
date ()
Sys.info ()[c("nodename", "user")]
commandArgs ()
rm (list = ls ())
R.version.string ## "R version 3.3.2 (2016-10-31)"

# libraries
library(dplyr)

# set the working directory
setwd("~/Documents/STRs/VALIDATION/raw_data/")

# Load validation golden table data
val_data = read.csv("../STRVALIDATION_ALLDATA_2019-10-7_ALL_kibanez.tsv",
                    header = T,
                    sep = "\t",
                    stringsAsFactors = F)
dim(val_data)
# 639 18

# Load Pilot merged table
merged_maxCI_table_pilot = read.csv("./pilot_validation/merged/merged_validation_pilot_maxCI.tsv",
                              sep = "\t",
                              header = T,
                              stringsAsFactors = F)
dim(merged_maxCI_table_pilot)
# 294 11

merged_avg_table_pilot = read.csv("./pilot_validation/merged/merged_validation_pilot_avg.tsv",
                                  sep = "\t",
                                  header = T,
                                  stringsAsFactors = F)
dim(merged_avg_table_pilot)
# 296  11

merged_avg_table_pilot = merged_avg_table_pilot %>%
  select(gene, allele, list_samples)
colnames(merged_avg_table_pilot) = c("gene", "avg_repeat", "list_samples")
merged_avg_table_pilot = merged_avg_table_pilot %>%
  mutate(unique_id = paste(gene, avg_repeat, sep = "_"))
dim(merged_avg_table_pilot)
# 296  4

merged_maxCI_table_pilot = merged_maxCI_table_pilot %>%
  select(gene, allele, list_samples)
colnames(merged_maxCI_table_pilot) = c("gene", "maxCI_repeat", "list_samples")
merged_maxCI_table_pilot = merged_maxCI_table_pilot %>%
  mutate(unique_id = paste(gene, maxCI_repeat, sep = "_"))
dim(merged_maxCI_table_pilot)
# 294  4

merged_all_pilot = full_join(merged_avg_table_pilot,
                         merged_maxCI_table_pilot,
                         by = c("unique_id"))
dim(merged_all_pilot)
# 320  7

colnames(merged_all_pilot) = c("gene_avg", "avg_repeat", "list_samples_avg", "unique_id", "gene_max", "maxCI_repeat", "list_samples_maxCI")


# Load Research merged table
merged_maxCI_table_research = read.csv("../raw_data/research_validation/merged/merged_validation_research_maxCI.tsv",
                                       sep = "\t",
                                       header = T,
                                       stringsAsFactors = F)
dim(merged_maxCI_table_research)
# 869  11

merged_maxCI_table_research = merged_maxCI_table_research %>%
  select(gene, allele, list_samples)
colnames(merged_maxCI_table_research) = c("gene", "maxCI_repeat", "list_samples")
merged_maxCI_table_research = merged_maxCI_table_research %>%
  mutate(unique_id = paste(gene, maxCI_repeat, sep = "_"))
dim(merged_maxCI_table_research)
# 869  4

merged_avg_table_research = read.csv("../raw_data/research_validation/merged/merged_validation_research_avg.tsv",
                                     sep = "\t",
                                     header = T,
                                     stringsAsFactors = F)
dim(merged_avg_table_research)
# 861  11

merged_avg_table_research = merged_avg_table_research %>%
  select(gene, allele, list_samples)
colnames(merged_avg_table_research) = c("gene", "avg_repeat", "list_samples")
merged_avg_table_research = merged_avg_table_research %>%
  mutate(unique_id = paste(gene, avg_repeat, sep = "_"))
dim(merged_avg_table_research)
# 861  4

merged_all_research = full_join(merged_avg_table_research,
                                merged_maxCI_table_research,
                                by = c("unique_id"))

dim(merged_all_research)
# 1391  7

colnames(merged_all_research) = c("gene_avg", "avg_repeat", "list_samples_avg", "unique_id", "gene_max", "maxCI_repeat", "list_samples_maxCI")

# Here should be WESSEX data (It's not research nor pilot)!!


# Let's enrich the validation golden table with the max CI value for each expansion
val_data2 = val_data
val_data2$EH_a1_avg = rep('NA', length(val_data$EH_a1))
val_data2$EH_a2_avg = rep('NA', length(val_data$EH_a2))
val_data2$EH_a1_maxCI = rep('NA', length(val_data$EH_a1))
val_data2$EH_a2_maxCI = rep('NA', length(val_data$EH_a2))

# We will go through all val_data rows, independently, one by one
# We will distinguish them by `gene` and `platekey`
for (i in 1:length(val_data$loci)){
  locus = val_data$locus_bioinfo[i]
  platekey = trimws(val_data$LP_Number[i])
  
  # avg values (to double check)
  row_avg_research = merged_all_research %>%
    filter(gene_avg %in% locus, grepl(platekey, list_samples_avg)) %>%
    select(avg_repeat) %>% pull() %>% as.character() 
  
  l_samples_avg = merged_all_research %>%
    filter(gene_avg %in% locus, grepl(platekey, list_samples_avg)) %>%
    select(list_samples_avg) %>% pull()
  
  if (length(l_samples_avg) > 0){
    if (grepl(paste(paste("EH_", platekey, sep = ""), "", sep = ".vcf_x2"), l_samples_avg)){
      row_avg_research = c(row_avg_research, row_avg_research)
    }
  }
  
  # maxCI values
  row_maxCI_research = merged_all_research %>%
    filter(gene_max %in% locus, grepl(platekey, list_samples_maxCI)) %>%
    select(maxCI_repeat) %>% pull() %>% as.character()
  
  l_samples_max = merged_all_research %>%
    filter(gene_max %in% locus, grepl(platekey, list_samples_maxCI)) %>%
    select(list_samples_maxCI) %>% pull()
  
  if (length(l_samples_max) > 0){
    if (grepl(paste(paste("EH_", platekey, sep = ""), "", sep = ".vcf_x2"), l_samples_max)){
      row_maxCI_research = c(row_maxCI_research, row_maxCI_research)
    }
  }
  
  if (length(row_avg_research) > 0){
    if (length(row_avg_research) < 2){
      val_data2$EH_a1_avg[i] = row_avg_research
      val_data2$EH_a1_maxCI[i] = row_maxCI_research
    }else{
      # Check whether we do have mod 2, otherwise do `unique`
      if (length(row_avg_research) %% 2 != 0){
        row_avg_research = unique(row_avg_research)
      }

      if (length(row_maxCI_research) %% 2 != 0){
        row_maxCI_research = unique(row_maxCI_research)
      }
      
      val_data2$EH_a1_avg[i] = row_avg_research[1]
      val_data2$EH_a2_avg[i] = row_avg_research[2]
      val_data2$EH_a1_maxCI[i] = row_maxCI_research[1]
      val_data2$EH_a2_maxCI[i] = row_maxCI_research[2]
    }
  }else{
    # IT's PILOT
    
    # avg values (to double check)
    row_avg_research = merged_all_pilot %>%
      filter(gene_avg %in% locus, grepl(platekey, list_samples_avg)) %>%
      select(avg_repeat) %>% pull() %>% as.character()
    
    l_samples_avg = merged_all_pilot %>%
      filter(gene_avg %in% locus, grepl(platekey, list_samples_avg)) %>%
      select(list_samples_avg) %>% pull()
    
    if (length(l_samples_avg) > 0){
      if (grepl(paste(paste("EH_", platekey, sep = ""), "", sep = ".vcf_x2"), l_samples_avg)){
        row_avg_research = c(row_avg_research, row_avg_research)
      }
    }
    
    # maxCI values
    row_maxCI_research = merged_all_pilot %>%
      filter(gene_max %in% locus, grepl(platekey, list_samples_maxCI)) %>%
      select(maxCI_repeat)  %>% pull() %>% as.character()
    
    l_samples_max = merged_all_pilot %>%
      filter(gene_max %in% locus, grepl(platekey, list_samples_maxCI)) %>%
      select(list_samples_maxCI) %>% pull()
    
    if (length(l_samples_max) > 0){
      if (grepl(paste(paste("EH_", platekey, sep = ""), "", sep = ".vcf_x2"), l_samples_max)){
        row_maxCI_research = c(row_maxCI_research, row_maxCI_research)
      }
    }
    
    if (length(row_avg_research) > 0){
      if (length(row_avg_research) < 2){
        val_data2$EH_a1_avg[i] = row_avg_research
        val_data2$EH_a1_maxCI[i] = row_maxCI_research
      }else{
        if (length(row_avg_research) %% 2 != 0){
          row_avg_research = unique(row_avg_research)
        }
        
        if (length(row_maxCI_research) %% 2 != 0){
          row_maxCI_research = unique(row_maxCI_research)
        }
        
        val_data2$EH_a1_avg[i] = row_avg_research[1]
        val_data2$EH_a2_avg[i] = row_avg_research[2]
        val_data2$EH_a1_maxCI[i] = row_maxCI_research[1]
        val_data2$EH_a2_maxCI[i] = row_maxCI_research[2]
      }
    }
  }
}

# Write results into file
write.table(val_data2, "../../ANALYSIS/EHv2_avg_VS_EHv2_maxCI/STRVALIDATION_ALLDATA_2019-10-7_ALL_kibanez_EHv255_avg_VS_EHv255_maxCI.tsv", quote = F, row.names = F, col.names = T, sep = "\t")
