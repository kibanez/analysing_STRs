# Objective: Complete our validation table together with the max CI value estimated by EH for each allele
# EHv3.1.2
date ()
Sys.info ()[c("nodename", "user")]
commandArgs ()
rm (list = ls ())
R.version.string ## "R version 3.6.1 (2019-07-05)"

# libraries
library(dplyr)

# set the working directory
setwd("~/Documents/STRs/VALIDATION/raw_data/EH-v3.1.2/")

# Load validation golden table data
val_data = read.csv("../../STRVALIDATION_ALLDATA_2019-10-7_ALL_kibanez_enriched.tsv",
                    header = T,
                    sep = "\t",
                    stringsAsFactors = F)
dim(val_data)
# 638 20

# Load Pilot merged table
merged_maxCI_table_pilot = read.csv("pilot_validation/merged/merged_validation_pilot_maxCI_EHv3.1.2.tsv",
                              sep = "\t",
                              header = T,
                              stringsAsFactors = F)
dim(merged_maxCI_table_pilot)
# 584 12

merged_avg_table_pilot = read.csv("pilot_validation/merged/merged_validation_pilot_avg_EHv3.1.2.tsv",
                                  sep = "\t",
                                  header = T,
                                  stringsAsFactors = F)
dim(merged_avg_table_pilot)
# 564  12

merged_avg_table_pilot = merged_avg_table_pilot %>%
  select(gene, allele, list_samples)
colnames(merged_avg_table_pilot) = c("gene", "avg_repeat", "list_samples")

merged_avg_table_pilot = merged_avg_table_pilot %>%
  mutate(unique_id = paste(gene, avg_repeat, sep = "_"))
dim(merged_avg_table_pilot)
# 564  4

merged_maxCI_table_pilot = merged_maxCI_table_pilot %>%
  select(gene, allele, list_samples)
colnames(merged_maxCI_table_pilot) = c("gene", "maxCI_repeat", "list_samples")
merged_maxCI_table_pilot = merged_maxCI_table_pilot %>%
  mutate(unique_id = paste(gene, maxCI_repeat, sep = "_"))
dim(merged_maxCI_table_pilot)
# 584  4

merged_all_pilot = full_join(merged_avg_table_pilot,
                         merged_maxCI_table_pilot,
                         by = c("unique_id"))
dim(merged_all_pilot)
# 645  7

colnames(merged_avg_table_pilot) = c("gene_avg", "avg_repeat", "list_samples_avg", "unique_id")
colnames(merged_maxCI_table_pilot) = c("gene_max", "maxCI_repeat", "list_samples_maxCI", "unique_id")
colnames(merged_all_pilot) = c("gene_avg", "avg_repeat", "list_samples_avg", "unique_id", "gene_max", "maxCI_repeat", "list_samples_maxCI")


# Load Research merged table
merged_maxCI_table_research = read.csv("research_validation/merged/merged_validation_research_maxCI_EHv3.1.2.tsv",
                                       sep = "\t",
                                       header = T,
                                       stringsAsFactors = F)
dim(merged_maxCI_table_research)
# 1199  12

merged_maxCI_table_research = merged_maxCI_table_research %>%
  select(gene, allele, list_samples)
colnames(merged_maxCI_table_research) = c("gene", "maxCI_repeat", "list_samples")
merged_maxCI_table_research = merged_maxCI_table_research %>%
  mutate(unique_id = paste(gene, maxCI_repeat, sep = "_"))
dim(merged_maxCI_table_research)
# 1199  4

merged_avg_table_research = read.csv("research_validation/merged/merged_validation_research_avg_EHv3.1.2.tsv",
                                     sep = "\t",
                                     header = T,
                                     stringsAsFactors = F)
dim(merged_avg_table_research)
# 1178  12

merged_avg_table_research = merged_avg_table_research %>%
  select(gene, allele, list_samples)
colnames(merged_avg_table_research) = c("gene", "avg_repeat", "list_samples")
merged_avg_table_research = merged_avg_table_research %>%
  mutate(unique_id = paste(gene, avg_repeat, sep = "_"))
dim(merged_avg_table_research)
# 1178  4

merged_all_research = full_join(merged_avg_table_research,
                                merged_maxCI_table_research,
                                by = c("unique_id"))

dim(merged_all_research)
# 1905  7

colnames(merged_avg_table_research) = c("gene_avg", "avg_repeat", "list_samples_avg", "unique_id")
colnames(merged_maxCI_table_research) = c("gene_max", "maxCI_repeat", "list_samples_maxCI", "unique_id")
colnames(merged_all_research) = c("gene_avg", "avg_repeat", "list_samples_avg", "unique_id", "gene_max", "maxCI_repeat", "list_samples_maxCI")

# Here should be WESSEX data (It's not research nor pilot)!!
# Load Wessex data
merged_maxCI_table_wessex = read.csv("./wessex_validation/merged/merged_validation_wessex_maxCI_EHv3.1.2.tsv",
                                    sep = "\t",
                                    header = T,
                                    stringsAsFactors = F)
dim(merged_maxCI_table_wessex)
# 394 12

merged_avg_table_wessex = read.csv("./wessex_validation//merged/merged_validation_wessex_avg_EHv3.1.2.tsv",
                                  sep = "\t",
                                  header = T,
                                  stringsAsFactors = F)
dim(merged_avg_table_wessex)
# 393  12

merged_avg_table_wessex = merged_avg_table_wessex %>%
  select(gene, allele, list_samples)
colnames(merged_avg_table_wessex) = c("gene", "avg_repeat", "list_samples")
merged_avg_table_wessex = merged_avg_table_wessex %>%
  mutate(unique_id = paste(gene, avg_repeat, sep = "_"))
dim(merged_avg_table_wessex)
# 393  4

merged_maxCI_table_wessex = merged_maxCI_table_wessex %>%
  select(gene, allele, list_samples)
colnames(merged_maxCI_table_wessex) = c("gene", "maxCI_repeat", "list_samples")
merged_maxCI_table_wessex = merged_maxCI_table_wessex %>%
  mutate(unique_id = paste(gene, maxCI_repeat, sep = "_"))
dim(merged_maxCI_table_wessex)
# 394  4

merged_all_wessex = full_join(merged_avg_table_wessex,
                             merged_maxCI_table_wessex,
                             by = c("unique_id"))
dim(merged_all_wessex)
# 638  7

colnames(merged_avg_table_wessex) = c("gene_avg", "avg_repeat", "list_samples_avg", "unique_id")
colnames(merged_maxCI_table_wessex) = c("gene_max", "maxCI_repeat", "list_samples_maxCI", "unique_id")
colnames(merged_all_wessex) = c("gene_avg", "avg_repeat", "list_samples_avg", "unique_id", "gene_max", "maxCI_repeat", "list_samples_maxCI")

# Let's enrich the validation golden table with the max CI value for each expansion
val_data2 = val_data
val_data2$EH_a1_avg = rep('NA', length(val_data$EH_a1))
val_data2$EH_a2_avg = rep('NA', length(val_data$EH_a2))
val_data2$EH_a1_maxCI = rep('NA', length(val_data$EH_a1))
val_data2$EH_a2_maxCI = rep('NA', length(val_data$EH_a2))
val_data2$gender = rep('NA', length(val_data$EH_a2))

# Research environment clinical data
clinical_data = read.csv("/Users/kibanez/Documents/STRs/clinical_data/clinical_data/rd_genomes_all_data_250919.tsv",
                         sep = "\t",
                         header = T,
                         stringsAsFactors = F)
dim(clinical_data)
# 1056568  26

# Let's retrieve only platekey and gender
clinical_data = clinical_data %>% 
  filter(plate_key.x %in% unique(val_data$LP_Number)) %>%
  select(plate_key.x, participant_phenotypic_sex)

dim(clinical_data)
# 8376  2

# Pilot clinical data
pilot_clinical_data = read.csv("/Users/kibanez/Documents/STRs/clinical_data/pilot_clinical_data/pilot_cohort_clinical_data_4833_genomes_removingPanels_280919.tsv",
                               sep = "\t",
                               header = T,
                               stringsAsFactors = F)
dim(pilot_clinical_data)
# 4974  10

pilot_clinical_data = pilot_clinical_data %>%
  filter(plateKey %in% unique(val_data$LP_Number)) %>%
  select(plateKey, sex)
dim(pilot_clinical_data)
# 69  2

# We will go through all val_data rows, independently, one by one
# We will distinguish them by `gene` and `platekey`
for (i in 1:length(val_data$loci)){
  locus = trimws(val_data$loci[i])
  platekey = trimws(val_data$LP_Number[i])
  
  # avg values (to double check)
  row_avg_research = merged_avg_table_research %>%
    filter(gene_avg %in% toupper(locus), grepl(platekey, list_samples_avg)) %>%
    select(avg_repeat) %>% pull() %>% as.character() 
  
  l_samples_avg = merged_avg_table_research %>%
    filter(gene_avg %in% toupper(locus), grepl(platekey, list_samples_avg)) %>%
    select(list_samples_avg) %>% pull()
  
  if (length(l_samples_avg) > 0){
    if (any(grepl(paste(paste("EH_", platekey, sep = ""), "", sep = ".vcf_x2"), l_samples_avg))){
      row_avg_research = c(row_avg_research, row_avg_research)
    }
  }
  
  # maxCI values
  row_maxCI_research = merged_maxCI_table_research %>%
    filter(gene_max %in% toupper(locus), grepl(platekey, list_samples_maxCI)) %>%
    select(maxCI_repeat) %>% pull() %>% as.character()
  
  l_samples_max = merged_maxCI_table_research %>%
    filter(gene_max %in% toupper(locus), grepl(platekey, list_samples_maxCI)) %>%
    select(list_samples_maxCI) %>% pull()
  
  if (length(l_samples_max) > 0){
    if (any(grepl(paste(paste("EH_", platekey, sep = ""), "", sep = ".vcf_x2"), l_samples_max))){
      row_maxCI_research = c(row_maxCI_research, row_maxCI_research)
    }
  }
  
  if (length(row_avg_research) > 0){
    if (length(row_avg_research) < 2){
      val_data2$EH_a1_avg[i] = row_avg_research
      val_data2$EH_a1_maxCI[i] = row_maxCI_research
    }else{
      val_data2$EH_a1_avg[i] = row_avg_research[1]
      val_data2$EH_a2_avg[i] = row_avg_research[2]
      val_data2$EH_a1_maxCI[i] = row_maxCI_research[1]
      val_data2$EH_a2_maxCI[i] = row_maxCI_research[2]
      
      if (platekey %in% unique(clinical_data$plate_key.x)){
        val_data2$gender[i] = clinical_data %>% 
          filter(plate_key.x %in% platekey) %>% 
          select(participant_phenotypic_sex) %>%
          pull() %>%
          unique() %>%
          as.character() 
      }else{
        val_data2$gender[i] = "NA"
      }
      
    }
  }else{
    # IT's PILOT
    
    # avg values (to double check)
    row_avg_research = merged_avg_table_pilot %>%
      filter(gene_avg %in% toupper(locus), grepl(platekey, list_samples_avg)) %>%
      select(avg_repeat) %>% pull() %>% as.character()
    
    l_samples_avg = merged_avg_table_pilot %>%
      filter(gene_avg %in% toupper(locus), grepl(platekey, list_samples_avg)) %>%
      select(list_samples_avg) %>% pull()
    
    if (length(l_samples_avg) > 0){
      if (any(grepl(paste(paste("EH_", platekey, sep = ""), "", sep = ".vcf_x2"), l_samples_avg))){
        row_avg_research = c(row_avg_research, row_avg_research)
      }
    }
    
    # maxCI values
    row_maxCI_research = merged_maxCI_table_pilot %>%
      filter(gene_max %in% toupper(locus), grepl(platekey, list_samples_maxCI)) %>%
      select(maxCI_repeat)  %>% pull() %>% as.character()
    
    l_samples_max = merged_maxCI_table_pilot %>%
      filter(gene_max %in% toupper(locus), grepl(platekey, list_samples_maxCI)) %>%
      select(list_samples_maxCI) %>% pull()
    
    if (length(l_samples_max) > 0){
      if (any(grepl(paste(paste("EH_", platekey, sep = ""), "", sep = ".vcf_x2"), l_samples_max))){
        row_maxCI_research = c(row_maxCI_research, row_maxCI_research)
      }
    }
    
    if (length(row_avg_research) > 0){
      if (length(row_avg_research) < 2){
        val_data2$EH_a1_avg[i] = row_avg_research
        val_data2$EH_a1_maxCI[i] = row_maxCI_research
      }else{
        val_data2$EH_a1_avg[i] = row_avg_research[1]
        val_data2$EH_a2_avg[i] = row_avg_research[2]
        val_data2$EH_a1_maxCI[i] = row_maxCI_research[1]
        val_data2$EH_a2_maxCI[i] = row_maxCI_research[2]
        
        if (platekey %in% unique(pilot_clinical_data$plateKey)){
          val_data2$gender[i] = pilot_clinical_data %>% 
            filter(plateKey %in% platekey) %>% 
            select(sex) %>%
            pull() %>%
            unique() %>%
            as.character() 
        }else{
          val_data2$gender[i] = "NA"
        }
        
      }
    }else{
      # It's WESSEX
      
      # avg values (to double check)
      row_avg_research = merged_avg_table_wessex %>%
        filter(gene_avg %in% toupper(locus), grepl(platekey, list_samples_avg)) %>%
        select(avg_repeat) %>% pull() %>% as.character()
      
      l_samples_avg = merged_avg_table_wessex %>%
        filter(gene_avg %in% toupper(locus), grepl(platekey, list_samples_avg)) %>%
        select(list_samples_avg) %>% pull()
      
      if (length(l_samples_avg) > 0){
        if (any(grepl(paste(paste("EH_", platekey, sep = ""), "", sep = ".vcf_x2"), l_samples_avg))){
          row_avg_research = c(row_avg_research, row_avg_research)
        }
      }
      
      # maxCI values
      row_maxCI_research = merged_maxCI_table_wessex %>%
        filter(gene_max %in% toupper(locus), grepl(platekey, list_samples_maxCI)) %>%
        select(maxCI_repeat)  %>% pull() %>% as.character()
      
      l_samples_max = merged_maxCI_table_wessex %>%
        filter(gene_max %in% toupper(locus), grepl(platekey, list_samples_maxCI)) %>%
        select(list_samples_maxCI) %>% pull()
      
      if (length(l_samples_max) > 0){
        if (any(grepl(paste(paste("EH_", platekey, sep = ""), "", sep = ".vcf_x2"), l_samples_max))){
          row_maxCI_research = c(row_maxCI_research, row_maxCI_research)
        }
      }
      
      if (length(row_avg_research) > 0){
        if (length(row_avg_research) < 2){
          val_data2$EH_a1_avg[i] = row_avg_research
          val_data2$EH_a1_maxCI[i] = row_maxCI_research
        }else{
          val_data2$EH_a1_avg[i] = row_avg_research[1]
          val_data2$EH_a2_avg[i] = row_avg_research[2]
          val_data2$EH_a1_maxCI[i] = row_maxCI_research[1]
          val_data2$EH_a2_maxCI[i] = row_maxCI_research[2]
          
          if (platekey %in% unique(clinical_data$plate_key.x)){
            val_data2$gender[i] = clinical_data %>% 
              filter(plate_key.x %in% platekey) %>% 
              select(participant_phenotypic_sex) %>%
              pull() %>%
              unique() %>%
              as.character() 
          }else{
            val_data2$gender[i] = "NA"
          }
          
        }
      }
    }
  }
}

# Select the columns we want to work with (not all, they already are in the validation golden table)
val_data2 = val_data2 %>%
  select(locus_bioinfo, LP_Number, gender, EH_a1, EH_a2, experimental_a1, experimental_a2, EH_a1_avg, EH_a2_avg, EH_a1_maxCI, EH_a2_maxCI, classification)

# To make things easier, we want to automatise the classification of each expansion call
table_threshold_normal = read.csv("/Users/kibanez/git/analysing_STRs/threshold_largest_normal_reported_research.txt",
                              stringsAsFactors = F,
                              header = T,
                              sep = "\t")

dim(table_threshold_normal)
# 32  2
colnames(table_threshold_normal) = c("locus_bioinfo", "threshold_normal")

table_threshold_pathogenic = read.csv("/Users/kibanez/git/analysing_STRs/threshold_smallest_pathogenic_reported_research.txt",
                                stringsAsFactors = F,
                                header = T,
                                sep = "\t")
dim(table_threshold_pathogenic)
# 32  2
colnames(table_threshold_pathogenic) = c("locus_bioinfo", "threshold_pathogenic")

# Let's include thresholds in the main table
val_data2 = full_join(val_data2,
                      table_threshold_normal,
                      by = "locus_bioinfo")

val_data2 = full_join(val_data2,
                      table_threshold_pathogenic,
                      by = "locus_bioinfo")

# Let's transform the data now
val_data2$EH_a1_avg = as.numeric(val_data2$EH_a1_avg)
val_data2$EH_a2_avg = as.numeric(val_data2$EH_a2_avg)
val_data2$EH_a1_maxCI = as.numeric(val_data2$EH_a1_maxCI)
val_data2$EH_a2_maxCI = as.numeric(val_data2$EH_a2_maxCI)

# Ignore all NA's when making OR and AND operations
val_data2$EH_a1_avg[which(is.na(val_data2$EH_a1_avg))] = 0
val_data2$EH_a2_avg[which(is.na(val_data2$EH_a2_avg))] = 0
val_data2$EH_a1_maxCI[which(is.na(val_data2$EH_a1_maxCI))] = 0
val_data2$EH_a2_maxCI[which(is.na(val_data2$EH_a2_maxCI))] = 0


# Let's define NEW CLASSIFICATION FOR AVG VALUES
val_data2 = val_data2 %>%
  group_by(LP_Number, locus_bioinfo) %>%
  mutate(new_classification_avg = case_when(((EH_a1_avg > threshold_normal | EH_a2_avg > threshold_normal) & (experimental_a1 == "expanded" | experimental_a2 == "expanded")) ~ "TP",
                                            ((EH_a1_avg > threshold_normal | EH_a2_avg > threshold_normal) & (experimental_a1 == "normal" & experimental_a2 == "normal")) ~ "FP",
                                            ((EH_a1_avg < threshold_normal & EH_a2_avg < threshold_normal) & (experimental_a1 != "expanded" & experimental_a2 != "expanded")) ~ "TN",
                                            ((EH_a1_avg < threshold_normal & EH_a2_avg < threshold_normal) & (experimental_a1 == "expanded" | experimental_a2 == "expanded")) ~ "FN")) %>%
  as.data.frame()

# Let's define NEW CLASSIFICATION FOR MAX_CI VALUES
val_data2 = val_data2 %>%
  group_by(LP_Number, locus_bioinfo) %>%
  mutate(new_classification_maxCI = case_when(((EH_a1_maxCI > threshold_normal | EH_a2_maxCI > threshold_normal) & (experimental_a1 == "expanded" | experimental_a2 == "expanded")) ~ "TP",
                                            ((EH_a1_maxCI > threshold_normal | EH_a2_maxCI > threshold_normal) & (experimental_a1 != "expanded" & experimental_a2 != "expanded")) ~ "FP",
                                            ((EH_a1_maxCI < threshold_normal & EH_a2_maxCI < threshold_normal) & (experimental_a1 != "expanded" & experimental_a2 != "expanded")) ~ "TN",
                                            ((EH_a1_maxCI < threshold_normal & EH_a2_maxCI < threshold_normal) & (experimental_a1 == "expanded" | experimental_a2 == "expanded")) ~ "FN")) %>%
  as.data.frame()

# TODO we need to make a special thing for FXN (or future biallelic or recessive loci) 
# I'll leave this to do post creating the excel file, manually, since there are ~10 validations that are not correctly created...
# Write results into file
write.table(val_data2, "/Users/kibanez/Documents/STRs/ANALYSIS/pipeline_performance/EHv3_avg_VS_EHv3_maxCI/STRVALIDATION_ALLDATA_2019-10-7_ALL_kibanez_EHv312_avg_VS_EHv312_maxCI_checkFXN.tsv", 
            quote = F, 
            row.names = F, 
            col.names = T, 
            sep = "\t")
