# Objective: define algorithmia behind annotated case-control TSV file when selecting repeat-motifs that are enriched in cases rather than in controls
date ()
Sys.info ()[c("nodename", "user")]
commandArgs ()
rm (list = ls ())
R.version.string ## "R version 3.6.3 (2020-02-29)"

# libraries
library(dplyr)

# defining working directory
setwd("~/Documents/STRs/ANALYSIS/EHdn/EHdn-v0.9.0/case-control/analysis/")

# Loading data
cc_data = read.csv("./test_ALS/output/casecontrol_test_ALS_locus-based_annotated.tsv",
                   stringsAsFactors = F,
                   header = T,
                   sep = "\t")
dim(cc_data)
# 8832  9

# List of cases
l_cases = read.table("./test_ALS/input/main_6_cases.txt",
                     stringsAsFactors = F)
l_cases = l_cases$V1
length(l_cases)
# 6

# List of controls
l_controls = read.table("./test_ALS/input/main_11355_controls.txt", stringsAsFactors = F)
l_controls = l_controls$V1
length(l_controls)
# 11355

# Filtering criteria: pvalue <= 0.05 
cc_data_filtered = cc_data %>%
  filter(pvalue <= 0.05)
dim(cc_data_filtered)
# 186  9

# Filtering those genes or repeat-motifs that are significantly enriched in cases
# We need to analyse `counts` field
# Create cases and controls columns, and put there the % of how many are there??

# Define first columns

cc_data_filtered$list_cases = rep(".", length(cc_data_filtered$contig))
cc_data_filtered$list_controls = rep(".", length(cc_data_filtered$contig))
cc_data_filtered$perc_cases = rep(".", length(cc_data_filtered$contig))
cc_data_filtered$perc_controls = rep(".", length(cc_data_filtered$contig))

for (i in 1:length(cc_data_filtered$counts)){
  aver = strsplit(cc_data_filtered$counts[i], ",")[[1]]
  total_counts = length(aver)
  cases_l = c()
  controls_l = c()
  
  for (j in 1:total_counts){
    item = strsplit(aver[j], ":")[[1]][1]
    if (item %in% l_cases){
      cases_l = c(cases_l, item)
    }else if (item %in% l_controls){
      controls_l = c(controls_l, item)
    }
  }
  
  # Raw list of cases and controls
  if (length(cases_l) == 0){
    cc_data_filtered$list_cases[i] = '.'
  }else{
    cc_data_filtered$list_cases[i] = paste(unlist(as.character(cases_l)), collapse = ';')  
  }
  
  cc_data_filtered$list_controls[i] = paste(unlist(as.character(controls_l)), collapse = ';')
  
  # % of cases and controls
  if (length(cases_l) == 0){
    perc_cases = 0
  }else{
    perc_cases = length(cases_l) / length(l_cases)  
  }
  
  perc_controls = length(controls_l) / length(l_controls)
  
  #print(perc_cases)
  #print(perc_controls)
  
  cc_data_filtered$perc_cases[i] = perc_cases
  cc_data_filtered$perc_controls[i] = perc_controls
  
}

dim(cc_data_filtered)
# 186  13

# Remove the `counts` column
cc_data_filtered = cc_data_filtered[-9]

write.table(cc_data_filtered,
            "./test_ALS/output/casecontrol_test_ALS_locus-based_annotated_filtered_with_lists.tsv",
            sep = "\t",
            quote = F,
            row.names = F,
            col.names = T)

write.table(cc_data_filtered %>% select(contig, start, end, motif, gene, region, pvalue, bonf_pvalue, perc_cases, perc_controls),
            "./test_ALS/output/casecontrol_test_ALS_locus-based_annotated_filtered_removing_lists.tsv",
            sep = "\t",
            quote = F,
            row.names = F,
            col.names = T)
