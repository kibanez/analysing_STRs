# Objective: compute and analyse how many genomes present an expansion beyond premut and patho thrsehold
# After visual inspection
date ()
Sys.info ()[c("nodename", "user")]
commandArgs ()
rm (list = ls ())
R.version.string ## "R version 3.6.3

# libraries
library(dplyr); packageDescription ("dplyr", fields = "Version") #"0.8.3"

# set working dire
setwd("~/Documents/STRs/ANALYSIS/population_research/PAPER/expanded_genomes_main_pilot/feb2021/")

# Load pathogenic table
patho_table = read.csv("./beyond_full-mutation/13_loci_beyond__pathogenic_cutoff_38_EHv322_92K_population_15M.tsv",
                       stringsAsFactors = F,
                       header = T,
                       sep = "\t")

# Compute how many expansions has each platekey
patho_table = patho_table %>%
  group_by(platekey, Final.decision) %>%
  mutate(number_exp = n()) %>%
  ungroup() %>%
  as.data.frame()

# See which ones are `number_exp` == 2 and `Final.decision` = Yes
patho_table %>% filter(number_exp == 2, Final.decision %in% "Yes")

# Load premutation table (FMR1 not yet)
premut_table = read.csv("./beyond_premut/13loci_beyond_premut_cutoff_to_review_VGD_enriched_pathoFinalDecision_100621.tsv",
                        stringsAsFactors = F,
                        header = T,
                        sep = "\t")

# Compute how many expansions has each platekey
premut_table  = premut_table %>%
  group_by(platekey, Final.decision) %>%
  mutate(number_exp = n()) %>%
  ungroup() %>%
  as.data.frame()

# See which ones are `number_exp` == 2 and `Final.decision` = Yes
premut_table %>% filter(number_exp == 2, Final.decision %in% "Yes")

