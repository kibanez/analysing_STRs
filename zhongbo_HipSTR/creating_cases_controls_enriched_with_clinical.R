# Objective: create a table of cases and control genomes (cases = neuro-like, controls = not neuro) enriched with clinical data across HipSTR genes
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
setwd("~/Documents/STRs/data/research/batch_december2020/output_EHv3.2.2_vcfs/merged/")

# load clinical data
clin_data = read.csv("~/Documents/STRs/clinical_data/clinical_data/Main_RE_V12_and_Pilot_programmes.tsv",
                     sep = "\t",
                     stringsAsFactors = F)
dim(clin_data)
# 2472865  26

# Create a dataframe, with `platekey` and `type` being: case, control or pseudocontrol
# case -> proband RD and recruited under Neurological
# control -> proband RD not neuro OR cancer, year of birth >40
# pseudocase -> not proband, RD not affected but recruited in a family under neuro
# pseudocontrol -> not proband, RD not affected but recruited in a family under NOT neuro, year of birth >40

l_families = clin_data %>%
  filter(grepl("Neuro", diseasegroup_list, ignore.case = T)) %>%
  select(rare_diseases_family_id) %>%
  unique() %>%
  pull()
length(l_families)
# 14717

l_cases = clin_data %>%
  filter(rare_diseases_family_id %in% l_families, participant_type %in% "Proband") %>%
  select(platekey) %>%
  unique() %>%
  pull()
length(l_cases)
# 14814

# QC
clin_data %>% filter(platekey %in% l_cases) %>% select(biological_relationship_to_proband)%>% table()
#Maternal Cousin Sister                 Mother                    N/A                Proband 
#30                     96                1358777                    582 

l_maternal_cousin_sister = clin_data %>% filter(platekey %in% l_cases, biological_relationship_to_proband %in% "Maternal Cousin Sister") %>% select(platekey) %>% unique() %>% pull()
l_Mother = clin_data %>% filter(platekey %in% l_cases, biological_relationship_to_proband %in% "Mother") %>% select(platekey) %>% unique() %>% pull()
l_cases = l_cases[-which(l_cases %in% c(l_maternal_cousin_sister, l_Mother))]
clin_data %>% filter(platekey %in% l_cases) %>% select(biological_relationship_to_proband)%>% table()
#N/A Proband 
#1358761     582

l_controls = clin_data %>%
  filter((programme %in% "Cancer" & year_of_birth <= "1981") | 
           ((!rare_diseases_family_id %in% l_families) & participant_type %in% "Proband" & year_of_birth <= "1981")) %>%
  select(platekey) %>%
  unique() %>%
  pull()
length(l_controls)
# 39822

clin_data %>% filter(platekey %in% l_controls) %>% select(biological_relationship_to_proband)%>% table()
#N/A Proband 
#185326     885 

clin_data %>% filter(platekey %in% l_controls) %>% select(year_of_birth)%>% pull() %>% as.integer() %>% summary()
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#1917    1954    1964    1963    1973    2018

l_pseudocases = clin_data %>%
  filter(rare_diseases_family_id %in% l_families, 
         participant_type %in% "Relative") %>%
  select(platekey) %>%
  unique() %>%
  pull()
length(l_pseudocases)
# 19360

# QC
clin_data %>% filter(platekey %in% l_pseudocases) %>% select(biological_relationship_to_proband)%>% table()
# OK

l_pseudocontrols = clin_data %>%
  filter(!rare_diseases_family_id %in% l_families, 
         participant_type %in% "Relative", 
         year_of_birth <= "1981") %>%
  select(platekey) %>%
  unique() %>%
  pull()
length(l_pseudocontrols)
# 16216

clin_data %>% filter(platekey %in% l_pseudocontrols) %>% select(biological_relationship_to_proband)%>% table()

# All 3 lists should be different
length(intersect(l_cases, l_controls))
# 0
length(intersect(l_pseudocontrols, l_controls))
# 0
length(intersect(l_pseudocases, l_cases))
# 0
length(intersect(l_pseudocases, l_controls))
# 0 

df_cases = data.frame(platekey = l_cases, type = rep("case", length(l_cases)))
df_controls = data.frame(platekey = l_controls, type = rep("control", length(l_controls)))
df_pseudocases = data.frame(platekey = l_pseudocases, type = rep("pseudocase", length(l_pseudocases)))
df_pseudocontrols = data.frame(platekey = l_pseudocontrols, type = rep("pseudocontrol", length(l_pseudocontrols)))

df_all = rbind(df_cases,
               df_controls,
               df_pseudocases,
               df_pseudocontrols)
df_all = unique(df_all)
dim(df_all)
# 90210  2

write.table(df_all, 
            "./analysis/table_cases_controls_90210_genomes_cases_controls_pseudoca_pseudoco.csv",
            quote = F,
            col.names = T,
            row.names = F,
            sep = ",")

# After running python generating_case_controls_from_VCF_EHv3.2.2.py, let's merge with this table
#--s
#/Volumes/KIKU_STRs/data/research/batch_december2020/output_EHv3.2.2_vcfs/
#  --specs
#/Users/kibanez/Documents/STRs/data/research/batch_december2020/output_EHv3.2.2_vcfs/variant_catalog_GRCh38_Zhongbo_december2020.json
#--o
#merged_93430_genomes_EHv322_batch_HipSTR_Zhongbo.tsv
#--O
#/Users/kibanez/Documents/STRs/data/research/batch_december2020/output_EHv3.2.2_vcfs/merged/
merged_table = read.csv("~/Documents/STRs/data/research/batch_december2020/output_EHv3.2.2_vcfs/merged/merged_93430_genomes_EHv322_batch_HipSTR_Zhongbo.tsv",
                        stringsAsFactors = F, 
                        sep = "\t")
dim(merged_table)
# 18176388  5
colnames(merged_table) = c("platekey", "gene", "a1", "a2", "coverage")

l_platekeys = unique(merged_table$platekey)
length(l_platekeys)
# 93430
l_genes = unique(merged_table$gene)
length(l_genes)
# 197

#Create `min_allele` and `max_allele`
merged_table = merged_table %>%
  group_by(platekey, gene) %>%
  mutate(min_allele = min(a1,a2),
         max_allele = max(a1,a2)) %>%
  ungroup() %>%
  as.data.frame()

merged_table = merged_table %>%
  select(platekey, gene, min_allele, max_allele)

cc_table = data.frame()
for(i in 1:length(l_platekeys)){
  df_aux = merged_table %>% 
    filter(platekey %in% l_platekeys[i]) 
  platekey_type = df_all %>% filter(platekey %in% l_platekeys[i]) %>% select(type) %>% pull() %>% as.character()
  if (length(platekey_type) != 1){
    platekey_type = "NA"
  }
  
  itziar = pivot_wider(df_aux, names_from = gene, values_from = c(min_allele, max_allele)) %>% as.data.frame()
  itziar$type = platekey_type
  
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
# 93430  396

write.table(cc_table, "./cc_table_93430_genomes_197_genes.tsv", sep = "\t", quote = F, row.names = F, col.names = T)

# The big data is too big to enrich clinically
# Let's do it by splitting it in smaller chunks
cc_table_main = read.csv("./cc_table_93430_genomes_197_genes.tsv",
                    sep = "\t",
                    stringsAsFactors = F)
dim(cc_table_main)
# 93430  396

to_enrich = clin_data %>%
  select(platekey, rare_diseases_family_id, participant_phenotypic_sex, year_of_birth, programme, diseasegroup_list, diseasesubgroup_list, hpo_list)

# Enrich it with gender, age, onset, disease_group, diseaes_subgroup, programme, hpo_terms - by gene (otherwise R will abort)
l_genes = unique(merged_table$gene)
for (i in 1:length(l_genes)){
  min_allele_column = paste("min_allele", l_genes[i], sep = "_")
  max_allele_column = paste("max_allele", l_genes[i], sep = "_")
  
  cc_table = cc_table_main %>%
    select(platekey, type, all_of(min_allele_column), all_of(max_allele_column))
  
  cc_table = left_join(cc_table,
                       to_enrich,
                       by = "platekey")
  cc_table = unique(cc_table)
  dim(cc_table)
  
  # which duplicated?
  platekeys_duplicated = which(duplicated(cc_table$platekey))
  l_platekeys_duplicated = cc_table$platekey[platekeys_duplicated]
  
  # They all have `year_of_birth` with different values...as we are going to use cut-off for age, let's keep the largest one (younger if they are not)
  year_messed_up = cc_table %>%
    filter(platekey %in% l_platekeys_duplicated) %>%
    unique()
  dim(year_messed_up)
  # 52 35
  
  year_messed_up = year_messed_up %>%
    group_by(platekey) %>%
    mutate(youngest_year = max(year_of_birth)) %>%
    ungroup() %>%
    as.data.frame()
  
  year_messed_up = year_messed_up[,-7]
  year_messed_up = unique(year_messed_up)
  
  colnames(year_messed_up)[11] = "year_of_birth"
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

  output_file_gene = paste("table_cases_controls_93425_genomes_enriched_with_some_clinical_data", merged_table$gene[i], sep = "_")
  output_file_gene = paste(output_file_gene, ".tsv", sep = "")
  if (grepl("/", output_file_gene)){
    output_file_gene = gsub("/", "_", output_file_gene)
  }
  write.table(cc_table, output_file_gene, sep = "\t", quote = F, col.names = T, row.names = F)
}


