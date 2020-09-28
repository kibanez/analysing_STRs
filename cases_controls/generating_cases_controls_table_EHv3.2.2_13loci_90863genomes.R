# Objective: create a single table including all 13 loci, and 90,863 genomes for which we have repeat-size for all these genes
# Enrich with relatedness (`Yes` or `No`)
# Enrich with population (batch1 AND batch2)
# For EHv3.2.2
date ()
Sys.info ()[c("nodename", "user")]
commandArgs ()
rm (list = ls ())
R.version.string ## "R version 3.6.1 (2019-07-05)"

# libraries
library(dplyr)

# Set working directory
setwd("~/Documents/STRs/ANALYSIS/cases_controls/batch_march/EHv322/")

# Load data
merged_data = read.csv("/Users/kibanez/Documents/STRs/data/research/batch_march2020/output_EHv3.2.2/merged/merged_92663_genomes_EHv3.2.2.tsv",
                       sep = '\t',
                       header = T,
                       stringsAsFactors = F)
dim(merged_data)
# 8560  12

# Data from RE rather than from Catalog (this clinical data has been retrieved from RE on Sept 2019)
# Pilot
pilot_clin_data = read.csv("~/Documents/STRs/clinical_data/pilot_clinical_data/pilot_cohort_clinical_data_4833_genomes_removingPanels_280919.tsv",
                           sep = "\t",
                           stringsAsFactors = FALSE,
                           header = TRUE)
dim(pilot_clin_data)
# 4974  10

# Let´s put all panel names into 1 single string splitted by ','
list_panels_pilot = pilot_clin_data %>% group_by(gelID) %>% summarise(panel_list = toString(unique(specificDisease))) %>% ungroup() %>% as.data.frame()
dim(list_panels_pilot)
# 4833  2

pilot_clin_data = left_join(pilot_clin_data,
                            list_panels_pilot,
                            by = "gelID")
# Remove specificDisease
pilot_clin_data = pilot_clin_data[,-8]
pilot_clin_data = unique(pilot_clin_data)
dim(pilot_clin_data)
# 4833  10

# Let's enrich with popu data
pilot_popu_table = read.csv("~/Documents/STRs/ANALYSIS/population_research/PILOT_ANCESTRY/FINE_GRAINED_RF_classifications_incl_superPOP_prediction_final20191216.csv",
                            stringsAsFactors = F,
                            header = T,
                            sep = ",")
dim(pilot_popu_table)
# 4821  44

pilot_clin_data = left_join(pilot_clin_data,
                            pilot_popu_table %>% select(ID, bestGUESS_sub_pop, bestGUESS_super_pop, self_reported),
                            by = c("plateKey"="ID"))

pilot_clin_data = unique(pilot_clin_data)
dim(pilot_clin_data)
# 4834  13

# Enrich pilot clinical data with disease group and disease subgroup
pheno_table = read.csv("~/Documents/STRs/clinical_data/pilot_clinical_data/phenotyping_v140_2019-09-13_15-26-02.tsv",
                       stringsAsFactors = F,
                       header = T,
                       sep = "\t")
pheno_table = unique(pheno_table)
dim(pheno_table)
# 106  3

pilot_clin_data = left_join(pilot_clin_data,
                            pheno_table,
                            by = c("panel_list" = "specific_disease"))
pilot_clin_data = unique(pilot_clin_data)
dim(pilot_clin_data)
# 4834  15

# Main
clin_data = read.table("~/Documents/STRs/clinical_data/clinical_data/rd_genomes_all_data_080920_V10.tsv",
                       sep = "\t",
                       stringsAsFactors = FALSE, 
                       header = TRUE)
dim(clin_data)  
# 1184730  31

# Let´s put all panel names into 1 single string splitted by ','
list_panels = clin_data %>% group_by(participant_id) %>% summarise(panel_list = toString(unique(panel_name))) %>% ungroup() %>% as.data.frame()
dim(list_panels)
# 87028  2

# Let´s put all HPO terms into 1 single string splitted by ','
list_hpos = clin_data %>% group_by(participant_id) %>% summarise(hpo_list = toString(unique(hpo_term))) %>% ungroup() %>% as.data.frame()
dim(list_hpos)
# 87028  2

# Let's put all specific diseases into 1 single string splitted by ','
list_diseases = clin_data %>% group_by(participant_id) %>% summarise(diseases_list = toString(unique(normalised_specific_disease))) %>% ungroup() %>% as.data.frame()
dim(list_diseases)
# 87028  2

list_disease_group = clin_data %>% group_by(participant_id) %>% summarise(diseasegroup_list = toString(unique(disease_group))) %>% ungroup() %>% as.data.frame()
dim(list_disease_group)
# 87028  2

list_disease_subgroup = clin_data %>% group_by(participant_id) %>% summarise(diseasesubgroup_list = toString(unique(disease_sub_group))) %>% ungroup() %>% as.data.frame()
dim(list_disease_subgroup)
# 87028  2

# Remove the panels and hpo columns, and include the list of panels and hpo respectively
clin_data = clin_data %>% 
  select(participant_id, platekey, rare_diseases_family_id, biological_relationship_to_proband, genome_build, year_of_birth, participant_phenotypic_sex, programme, family_group_type, affection_status, best_guess_predicted_ancstry, self_reported)
dim(clin_data)
# 1184730  12

clin_data = left_join(clin_data,
                      list_diseases,
                      by = "participant_id")
dim(clin_data)
# 1184730  13

clin_data = left_join(clin_data,
                      list_disease_group,
                      by = "participant_id")
dim(clin_data)
# 1184730  14

clin_data = left_join(clin_data,
                      list_disease_subgroup,
                      by = "participant_id")
dim(clin_data)
# 1184730  15


clin_data = left_join(clin_data,
                      list_panels,
                      by = "participant_id")
dim(clin_data)
# 1184730  16

clin_data = left_join(clin_data,
                      list_hpos,
                      by = "participant_id")
dim(clin_data)
# 1184730   17

# Enrich clin_data with pilot_clin_data, keeping diff fields as `.`
colnames(pilot_clin_data) = c("participant_id", "platekey", "rare_diseases_family_id", "participant_phenotypic_sex", "biological_relationship_to_proband", "affection_status", "year_of_birth", "ageOfOnset", "qc_state", "diseases_list", "best_guess_predicted_ancstry", "bestGUESS_super_pop", "self_reported", "diseasesubgroup_list", "diseasegroup_list")

# Generate extra columns from clin data for pilot clin data
pilot_clin_data$genome_build = rep("GRCh37", length(pilot_clin_data$participant_id))
pilot_clin_data$programme = rep("RD Pilot", length(pilot_clin_data$participant_id))
pilot_clin_data$family_group_type = rep(".", length(pilot_clin_data$participant_id))
pilot_clin_data$panel_list = rep(".", length(pilot_clin_data$participant_id))
pilot_clin_data$hpo_list = rep(".", length(pilot_clin_data$participant_id))

# Select columns
pilot_clin_data = pilot_clin_data %>%
  select(participant_id, platekey, rare_diseases_family_id, biological_relationship_to_proband,
         genome_build, year_of_birth, participant_phenotypic_sex, programme,
         family_group_type, affection_status, best_guess_predicted_ancstry, self_reported,
         diseases_list, diseasegroup_list, diseasesubgroup_list, panel_list, hpo_list)

clin_data = rbind(clin_data,
                  pilot_clin_data)
dim(clin_data)
# 1189564  17

# As output
# We want to the following output - NOTE each row is a participant (!!!)
# Family ID
# Patient ID
# LP number
# Participant type (proband / mother / brother etc)
# HTT allele
# Disease Category
# Disease Subgroup
# Panels applied
# Sex
# Age at recruitment
# Age of onset
# no of family participants
# ethnic or popu
# unrelated (`yes` or `no`)

# Now in sep'20, we do have latest merged batch1 and batch2 data for population
popu_merged = read.csv("~/Documents/STRs/ANALYSIS/population_research/MAIN_ANCESTRY/popu_merged_batch1_batch2_79849_genomes.tsv",
                       sep = "\t",
                       stringsAsFactors = F,
                       header = F)
dim(popu_merged)
# 79849  2
colnames(popu_merged) = c("platekey", "superpopu")

l_unrelated_genomes= read.table("~/Documents/STRs/ANALYSIS/population_research/MAIN_ANCESTRY/batch2/l_unrelated_55603_genomes_batch2.txt",
                                stringsAsFactors = F)
l_unrelated_genomes = l_unrelated_genomes$V1
length(l_unrelated_genomes)
# 55603

# Let's enrich clin data with this
clin_data$unrelated = rep("No", length(clin_data$platekey))
clin_data$unrelated[which(clin_data$platekey %in% l_unrelated_genomes)] = "Yes"
dim(clin_data)
# 1189564 18

clin_data = unique(clin_data)
dim(clin_data)
# 111312  18

# Let's enrich with merged b1&b2 popu data
clin_data = left_join(clin_data,
                      popu_merged,
                      by ="platekey")
clin_data = unique(clin_data)
dim(clin_data)
# 111312  19

# clin_data contains population for merged (superpopu) and batch1 (best_guess_predicted_ancstry and self_reported) data
# let's enrich it with batch2 popu data
popu_batch2 = read.csv("~/Documents/STRs/ANALYSIS/population_research/MAIN_ANCESTRY/batch2/aggV2_M30K_60K_1KGP3_ancestry_assignment_probs_R9_08062020.tsv",
                       stringsAsFactors = F,
                       sep = " ",
                       header = T)
dim(popu_batch2)
# 78388  33

clin_data = left_join(clin_data,
                      popu_batch2 %>% select(plate_key, ancestry0_8),
                      by = c("platekey" = "plate_key"))
clin_data = unique(clin_data)
dim(clin_data)
# 111312  20

# We now want to have info for 13 loci
# and intersected 
l_genomes_across_selected_loci = read.table("~/Documents/STRs/ANALYSIS/population_research/PAPER/carriers/cc_pileup_100Kg/list_90863_unique_similar_genomes_across_13_loci.txt",
                            stringsAsFactors = F)
l_genomes_across_selected_loci = l_genomes_across_selected_loci$V1
length(l_genomes_across_selected_loci)
# 90863

clin_data_selected = clin_data %>%
  filter(platekey %in% l_genomes_across_selected_loci)
dim(clin_data_selected)
# 89821  19

#l_genes = sort(unique(merged_data$gene))
l_genes = c("AR", "ATN1", "ATXN1", "ATXN2", "ATXN3", "ATXN7", "CACNA1A", "C9ORF72", "DMPK", "FMR1", "FXN", "HTT", "TBP")

locus_data_genomic = c()
for (i in 1:length(l_genomes_across_selected_loci)){
  locus_data = merged_data %>% filter(grepl(l_genomes_across_selected_loci[i], list_samples),
                                      gene %in% l_genes) %>%
    select(allele, gene, list_samples)
  
  # We need to keep specifically with the sample_id, in order to know whether it's x2 or not
  for (num_allele in 1:length(locus_data$allele)){
    l_samples = strsplit(locus_data$list_samples[num_allele], ";")[[1]]
    index_sample = which(grepl(l_genomes_across_selected_loci[i], l_samples))
    if (grepl("x2", l_samples[index_sample])){
      locus_data = rbind(locus_data,
                         locus_data[num_allele,])
    }
  }# num_alelle
  
  # Remove `list_samples` from locus_data
  locus_data = locus_data %>% select(allele, gene)
  locus_data$platekey = rep(l_genomes_across_selected_loci[i], length(locus_data$allele))
  
  # Reformat df
  ar_alleles = locus_data %>% filter(gene %in% "AR") %>% t() %>% as.data.frame() 
  # chrX: if there is only one allele
  if (dim(ar_alleles)[2] < 2){
    ar_alleles = cbind(ar_alleles, rbind(NA, "AR", l_genomes_across_selected_loci[i]))
  }
  
  atn1_alleles = locus_data %>% filter(gene %in% "ATN1") %>% t() %>% as.data.frame() 
  atxn1_alleles = locus_data %>% filter(gene %in% "ATXN1") %>% t() %>% as.data.frame() 
  atxn2_alleles = locus_data %>% filter(gene %in% "ATXN2") %>% t() %>% as.data.frame() 
  atxn3_alleles = locus_data %>% filter(gene %in% "ATXN3") %>% t() %>% as.data.frame() 
  atxn7_alleles = locus_data %>% filter(gene %in% "ATXN7") %>% t() %>% as.data.frame() 
  cacna1a_alleles = locus_data %>% filter(gene %in% "CACNA1A") %>% t() %>% as.data.frame() 
  c9orf72_alleles = locus_data %>% filter(gene %in% "C9ORF72") %>% t() %>% as.data.frame() 
  dmpk_alleles = locus_data %>% filter(gene %in% "DMPK") %>% t() %>% as.data.frame() 
  fmr1_alleles = locus_data %>% filter(gene %in% "FMR1") %>% t() %>% as.data.frame() 
  # chrX: if there is only one allele
  if (dim(fmr1_alleles)[2] < 2){
    fmr1_alleles = cbind(fmr1_alleles, rbind(NA, "FMR1", l_genomes_across_selected_loci[i]))
  }
  
  fxn_alleles = locus_data %>% filter(gene %in% "FXN") %>% t() %>% as.data.frame() 
  htt_alleles = locus_data %>% filter(gene %in% "HTT") %>% t() %>% as.data.frame() 
  tbp_alleles = locus_data %>% filter(gene %in% "TBP") %>% t() %>% as.data.frame() 
  
  all_alleles = cbind(ar_alleles[1,], 
                      atn1_alleles[1,],
                      atxn1_alleles[1,],
                      atxn2_alleles[1,],
                      atxn3_alleles[1,],
                      atxn7_alleles[1,],
                      cacna1a_alleles[1,],
                      c9orf72_alleles[1,],
                      dmpk_alleles[1,],
                      fmr1_alleles[1,],
                      fxn_alleles[1,],
                      htt_alleles[1,],
                      tbp_alleles[1,])
  colnames(all_alleles) = c("AR_a1", "AR_a2",
                            "ATN1_a1", "ATN1_a2",
                            "ATXN1_a1", "ATXN1_a2",
                            "ATXN2_a1", "ATXN2_a2",
                            "ATXN3_a1", "ATXN3_a2",
                            "ATXN7_a1", "ATXN7_a2",
                            "CACNA1A_a1", "CACNA1A_a2",
                            "C9ORF72_a1", "C9ORF72_a2",
                            "DMPK_a1", "DMPK_a2",
                            "FMR1_a1", "FMR1_a2",
                            "FXN_a1", "FXN_a2",
                            "HTT_a1", "HTT_a2",
                            "TBP_a1", "TBP_a2")
  all_alleles$platekey = rep(l_genomes_across_selected_loci[i], length(all_alleles$AR_a1))
  
  locus_data_genomic = rbind(locus_data_genomic, all_alleles)
}

dena = locus_data_genomic

# Now that we have genomic (STR) data for this 90,863 genomes, let's merge them with clinical data
locus_data_new = left_join(locus_data_genomic,
                           clin_data_selected,
                           by = "platekey")
locus_data_new = unique(locus_data_new)
dim(locus_data_new)
# 91857  47

write.table(locus_data_new, "table_13_loci_across_90863_genomes_each_row_allele_EHv3.2.2_90863platekeys_88827pids.tsv", sep = "\t", quote = F, row.names = F, col.names = T)