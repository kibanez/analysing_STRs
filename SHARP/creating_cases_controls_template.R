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
# control -> RD not affected and not neuro
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
  filter(programme %in% "Cancer" | (!rare_diseases_family_id %in% l_families & affection_status %in% c("Unaffected", "NotAffected"))) %>%
  select(platekey) %>%
  unique() %>%
  pull()
length(l_controls)
# 50696

l_pseudocontrols = clin_data %>%
  filter(rare_diseases_family_id %in% l_families, affection_status %in% c("Unaffected", "NotAffected")) %>%
  select(platekey) %>%
  unique() %>%
  pull()
length(l_pseudocontrols)
# 16859

# All 3 lists should be different
length(intersect(l_cases, l_controls))
# 1 - "LP2000905-DNA_A01"
length(intersect(l_pseudocontrols, l_controls))
# 1 - "LP2000905-DNA_A01"
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
# 50695
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
# 83768  2

write.table(df_all, 
            "./analysis/table_cases_controls_83768_genomes.csv",
            quote = F,
            col.names = T,
            row.names = F,
            sep = ",")