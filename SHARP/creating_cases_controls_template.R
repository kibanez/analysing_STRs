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











