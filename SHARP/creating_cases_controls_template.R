# Objective: create a template of cases and control genomes (cases = neuro-like, controls = not neuro)
date ()
Sys.info ()[c("nodename", "user")]
commandArgs ()
rm (list = ls ())
R.version.string ## "R version 3.6.3 (2020-02-29)"

# libraries
library(ggplot2); packageDescription ("ggplot2", fields = "Version") #"3.2.1"
library(grid); packageDescription ("grid", fields = "Version") #"3.6.1"
library(gridExtra); packageDescription ("gridExtra", fields = "Version") #"2.3"
library(reshape2); packageDescription ("reshape2", fields = "Version") #"1.4.3"
require(dplyr); packageDescription ("dplyr", fields = "Version") #"0.8.3"
require(tidyr);packageDescription ("tidyr", fields = "Version") #"1.0.2"

# Set working dir
setwd("~/Documents/STRs/ANALYSIS/SHARP/EHdn_Parkinson/")

# load clinical data
clin_data = read.csv("~/Documents/STRs/clinical_data/clinical_data/Main_RE_V11_and_Pilot_programmes.tsv",
                     sep = "\t",
                     stringsAsFactors = F)
dim(clin_data)
# 2444984 24

# Create a dataframe, with `platekey` and `type` being: case, control or pseudocontrol
# case -> RD affected OR proband and recruited under Neurological 
# control -> RD not neuro (even if they are affected in other diseases)
# pseudocontrol -> RD not affected but recruited in a family under neuro

l_families = clin_data %>%
  filter(grepl("Neuro", diseasegroup_list, ignore.case = T)) %>%
  select(rare_diseases_family_id) %>%
  unique() %>%
  pull()
length(l_families)
# 14421

l_cases = clin_data %>%
  filter(rare_diseases_family_id %in% l_families, (biological_relationship_to_proband %in% "N/A" | affection_status %in% "Affected")) %>%
  select(platekey) %>%
  unique() %>%
  pull()
length(l_cases)
# 16224

l_controls = clin_data %>%
  filter(programme %in% "Cancer" | !rare_diseases_family_id %in% l_families) %>%
  select(platekey) %>%
  unique() %>%
  pull()
length(l_controls)
# 76364

l_pseudocontrols = clin_data %>%
  filter(rare_diseases_family_id %in% l_families, affection_status %in% c("Unaffected", "NotAffected")) %>%
  select(platekey) %>%
  unique() %>%
  pull()
length(l_pseudocontrols)
# 16859

# All 3 lists should be different
length(intersect(l_cases, l_controls))
# 2
length(intersect(l_pseudocontrols, l_controls))
# 2
length(intersect(l_cases, l_pseudocontrols))
# 5 - "LP2000905-DNA_A01" "LP3001339-DNA_E06" "LP3000032-DNA_D11" "LP3001180-DNA_C06" "LP3000448-DNA_E10"

# All these platekeys have 2 different and incongruent values for affection status: affected and unaffected - let's remove them from now
to_remove = c("LP2000905-DNA_A01","LP3001339-DNA_E06","LP3000032-DNA_D11","LP3001180-DNA_C06","LP3000448-DNA_E10")
index_to_remove_cases = pmatch(to_remove, l_cases)
index_to_remove_controls = pmatch(to_remove, l_controls)
index_to_remove_controls = index_to_remove_controls[!is.na(index_to_remove_controls)]
index_to_remove_pseudocontrols = pmatch(to_remove, l_pseudocontrols)

l_cases = l_cases[-index_to_remove_cases]
l_controls = l_controls[-index_to_remove_controls]
l_pseudocontrols = l_pseudocontrols[-index_to_remove_pseudocontrols]
length(l_cases)
# 16219
length(l_controls)
# 76362
length(l_pseudocontrols)
# 16854

df_cases = data.frame(platekey = l_cases, type = rep("case", length(l_cases)))
df_controls = data.frame(platekey = l_controls, type = rep("control", length(l_controls)))
df_pseudocontrols = data.frame(platekey = l_pseudocontrols, type = rep("pseudocontrol", length(l_pseudocontrols)))

df_all = rbind(df_cases,
               df_controls,
               df_pseudocontrols)
df_all = unique(df_all)
dim(df_all)
# 109435  2

write.table(df_all, 
            "./analysis/table_cases_controls_109435_genomes.csv",
            quote = F,
            col.names = T,
            row.names = F,
            sep = ",")

# After running python generating_case_controls_from_VCF_EHv3.2.2.py, let's merge with this table
# --s
#/Users/kibanez/Documents/STRs/ANALYSIS/SHARP/EHdn_Parkinson/output_EHv3.2.2_vcfs
#--specs
#/Users/kibanez/Documents/STRs/ANALYSIS/SHARP/EHdn_Parkinson/specs/EH_hg38_9LociForReplication_25March21_GRCh38.json
#--o
#merged_93425_genomes_EHv322_batch_Sharp_and_Parkinsonism.tsv
#--O
#/Users/kibanez/Documents/STRs/ANALYSIS/SHARP/EHdn_Parkinson/output_EHv3.2.2_vcfs/merged
merged_table = read.csv("~/Documents/STRs/ANALYSIS/SHARP/EHdn_Parkinson/output_EHv3.2.2_vcfs/merged/merged_93425_genomes_EHv322_batch_Sharp_and_Parkinsonism.tsv",
                        stringsAsFactors = F, 
                        sep = "\t")
dim(merged_table)
# 829327  5
colnames(merged_table) = c("platekey", "gene", "a1", "a2", "coverage")

l_platekeys = unique(merged_table$platekey)
length(l_platekeys)
# 93425
l_genes = unique(merged_table$gene)
length(l_genes)
# 9

cc_table = data.frame()
for(i in 1:length(l_platekeys)){
  df_aux = merged_table %>% filter(platekey %in% l_platekeys[i])
  platekey_type = df_all %>% filter(platekey %in% l_platekeys[i]) %>% select(type) %>% pull() %>% as.character()
  if (length(platekey_type != 1)){
    platekey_type = "NA"
  }
  
  itziar = pivot_wider(df_aux, names_from = gene, values_from = c(a1, a2, coverage)) %>% as.data.frame()

  # Check whether all genes are genotyped
  colnames_cc = colnames(cc_table)
  colnames_itziar = colnames(itziar)
  colnames_diff = setdiff(colnames_cc, colnames_itziar)
  if (length(colnames_diff) > 0){
    df_new_columns = data.frame()
    for(j in 1:length(colnames_diff)){
      df_new_columns = rbind(df_new_columns,
                             cbind(assign(paste0(colnames_diff[j], ''), colnames_diff[j]), "NA"))
    }
    df_new_columns = as.data.frame(t(df_new_columns))
    
    names(df_new_columns) = df_new_columns %>% slice(1) %>% unlist()
    df_new_columns <- df_new_columns %>% slice(-1)
    
    # Append to itziar
    itziar = cbind(itziar,
                   df_new_columns)
  }
  
  itziar = itziar[ , order(names(itziar))]
  cc_table = rbind(cc_table,
                   itziar)
}

dim(cc_table)
# 93425  28

write.table(cc_table, "./analysis/cc_table_93425_genomes_9_genes.tsv", sep = "\t", quote = F, row.names = F, col.names = T)

# Enrich it with gender, age, onset, disease_group, diseaes_subgroup, programme, hpo_terms
to_enrich = clin_data %>%
  select(platekey, participant_phenotypic_sex, year_of_birth, programme, diseasegroup_list, diseasesubgroup_list, hpo_list)

cc_table = left_join(cc_table,
                     to_enrich,
                     by = "platekey")
cc_table = unique(cc_table)
dim(cc_table)
# 93451 34

# which duplicated?
platekeys_duplicated = which(duplicated(cc_table$platekey))
length(platekeys_duplicated)
# 26
l_platekeys_duplicated = cc_table$platekey[platekeys_duplicated]

# They all have `year_of_birth` with different values...as we are going to use cut-off for age, let's keep the largest one (younger if they are not)
year_messed_up = cc_table %>%
  filter(platekey %in% l_platekeys_duplicated) %>%
  unique()
dim(year_messed_up)
# 52 34

year_messed_up = year_messed_up %>%
  group_by(platekey) %>%
  mutate(youngest_year = max(year_of_birth)) %>%
  ungroup() %>%
  as.data.frame()

year_messed_up = year_messed_up[,-30]
year_messed_up = unique(year_messed_up)
dim(year_messed_up)
# 30 34

colnames(year_messed_up)[34] = "year_of_birth"
year_messed_up = year_messed_up[colnames(cc_table)]  

# remove the duplicated ones
to_remove = year_messed_up$platekey[which(duplicated(year_messed_up$platekey))]

cc_table = cc_table %>%
  filter(!platekey %in% l_platekeys_duplicated)

cc_table = rbind(cc_table,
                 year_messed_up)

cc_table = cc_table %>%
  filter(!platekey %in% to_remove)
cc_table = unique(cc_table)
dim(cc_table)
# 93422 34

write.table(cc_table, "./analysis/table_cases_controls_93422_genomes_enriched_with_some_clinical_data.tsv", sep = "\t", quote = F, col.names = T, row.names = F)
