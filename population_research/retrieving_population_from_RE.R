# Objective: retrieve the ancestry estimation for each genome GEL research group team has estimated so far (with Plink)
date ()
Sys.info ()[c("nodename", "user")]
commandArgs ()
rm (list = ls ())
R.version.string ## "R version 3.3.2 (2016-10-31)"

# libraries
library(dplyr)

# set the working directory
setwd("~/Documents/STRs/clinical_data/clinical_data/")

popu_table = read.table("aggregate_gvcf_sample_stats_2019-10-03_22-07-31.tsv",
                        sep = "\t",
                        header = T,
                        stringsAsFactors = F)
dim(popu_table)
# 59356  51

# The research group at GEL has estimated with Plink the ancestry of 59,356 unrelated genomes (probands) RD and cancer germlines
popu_table = popu_table %>% select(participant_id, platekey, type, sample_type, participant_phenotypic_sex, pred_african_ancestries, pred_american_ancestries, pred_european_ancestries, pred_east_asian_ancestries, pred_south_asian_ancestries)
dim(popu_table)
# 59356  10

write.table(popu_table,
            "~/Documents/STRs/ANALYSIS/population_research/population_info_by_031019.tsv",
            sep = "\t",
            row.names = F,
            col.names = T,
            quote = F)

popu_table_enriched = popu_table %>% select(participant_id, platekey, type, sample_type, participant_phenotypic_sex, pc1, pc2, pc3, pc4, pc5, pc6, pc7, pc8, pc9, pc10, pred_african_ancestries, pred_american_ancestries, pred_european_ancestries, pred_east_asian_ancestries, pred_south_asian_ancestries)
dim(popu_table_enriched)
# 59356  20 

write.table(popu_table_enriched,
            "~/Documents/STRs/ANALYSIS/population_research/population_info_enriched_59356_by_031019.tsv",
            sep = "\t",
            row.names = F,
            col.names = T,
            quote = F)

