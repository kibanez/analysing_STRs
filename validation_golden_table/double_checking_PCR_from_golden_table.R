# Objective: we want to double check and be sure that the PCR results for each <PLATEKEY, LOCUS> are correct
# For that, we will compare the table from October 2019 (with AR PCR results updated from Liana/James) and compare them with the final table we do have for the papers 
date ()
Sys.info ()[c("nodename", "user")]
commandArgs ()
rm (list = ls ())
R.version.string ## "R version 3.6.1 (2019-07-05)"

# libraries
library(dplyr); packageDescription ("dplyr", fields = "Version") #"0.8.3"

# Set environment
setwd("/Users/kibanez/Documents/STRs/VALIDATION/")

# Load tables
paper_table = read.csv("EHv255_EHv312_validation_cohort_GEL.tsv",
                       stringsAsFactors = F,
                       header = T,
                       sep = "\t")
dim(paper_table)
# 638  28

# Note: this table contains 3 duplicates
orig_table = read.csv("STRVALIDATION_ALLDATA_2019-10-7_ALL_kibanez_enriched.tsv",
                      stringsAsFactors = F,
                      header = T,
                      sep = "\t")
dim(orig_table)
# 638  20

new_paper_table = paper_table
new_paper_table$exp_PCR_a1_checked = rep("", length(new_paper_table$locus))
new_paper_table$exp_PCR_a2_checked = rep("", length(new_paper_table$locus))

for(i in 1:length(new_paper_table$locus)){
  locus = new_paper_table$locus[i]
  platekey = new_paper_table$LP_number[i]
  
  new_exp_PCR_a1 = orig_table %>% filter(LP_Number %in% platekey, loci %in% locus) %>% select(STR_a1) %>% pull() 
  new_exp_PCR_a2 = orig_table %>% filter(LP_Number %in% platekey, loci %in% locus) %>% select(STR_a2) %>% pull() 
  
  if (length(new_exp_PCR_a1) > 1){
    new_exp_PCR_a1 = new_exp_PCR_a1[1]
  }else if (length(new_exp_PCR_a1 )){
    next
  }
  if (length(new_exp_PCR_a2) > 1){
    new_exp_PCR_a2 = new_exp_PCR_a2[1]
  }

    new_paper_table$exp_PCR_a1_checked = new_exp_PCR_a1
  new_paper_table$exp_PCR_a2_checked = new_exp_PCR_a2
}

write.table(new_paper_table,
            "EHv255_EHv312_validation_cohort_GEL_PCR_estimates_double_checked.tsv",
            quote = F,
            row.names = F,
            col.names = T,
            sep = "\t")
