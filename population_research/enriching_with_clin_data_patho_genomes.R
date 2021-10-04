# Script to enrich with clinical data genomes that have an expansion beyond the full-mutation cut-off
# and see whether they can be considered in our cohort
date ()
Sys.info ()[c("nodename", "user")]
commandArgs ()
rm (list = ls ())
R.version.string ## "R version 3.6.1 (2019-07-05)"

# libraries
library(dplyr)

# Set working directory
setwd("~/Documents/STRs/ANALYSIS/population_research/PAPER/expanded_genomes_main_pilot/feb2021/enriched_with_clin_data/")

# Load clinical data
clin_data = read.csv("~/Documents/STRs/clinical_data/clinical_data/Main_RE_V12_and_Pilot_programmes.tsv", stringsAsFactors = F, header = T, sep= "\t")
dim(clin_data)
# 2472865  26

# Load genomes to enrich
genomes_to_enrich = read.csv("./genomes_beyond_patho_with_clinical_data.tsv",
                             stringsAsFactors = F,
                             header = T,
                             sep = "\t")
dim(genomes_to_enrich)
# 108  12 

l_genes = unique(genomes_to_enrich$locus)
for(i in 1:length(l_genes)){
  to_write = left_join(genomes_to_enrich %>% filter(locus %in% l_genes[i]),
                       clin_data %>% select(participant_id, platekey, rare_diseases_family_id, biological_relationship_to_proband, year_of_birth, participant_phenotypic_sex, participant_type, affection_status, diseases_list, diseasegroup_list, diseasesubgroup_list, panel_list, hpo_list),
                       by = "participant_id")
  to_write = unique(to_write)
  output_file = paste(l_genes[i], "enriched_with_clinical_data.tsv", sep = "_")
  write.table(to_write, output_file, row.names = F, col.names = T, sep = "\t", quote = F)
}
