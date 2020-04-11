# Objective: create subpopuation labels INP file for fastPHASE input
date ()
Sys.info ()[c("nodename", "user")]
commandArgs ()
rm (list = ls ())
R.version.string ## "R version 3.6.1 (2019-07-05)"

require(dplyr); packageDescription ("dplyr", fields = "Version") #"0.8.3"

# Working directory
setwd("~/Documents/STRs/ANALYSIS/haplotyping/AR/")

# Create the dataframe we will use to name or label each population
superpopulations = c("AFR","AFR","AFR","AFR","AFR","AFR","AFR", 
                     "AMR", "AMR","AMR","AMR",
                     "EUR","EUR","EUR","EUR","EUR", 
                     "EAS","EAS","EAS","EAS","EAS",
                     "SAS","SAS","SAS","SAS","SAS")
sub_populations = c("ESN", "YRI", "GWD", "LWK", "MSL", "ACB", "ASW", 
                    "MXL", "PUR", "PEL", "CLM",
                    "IBS", "TSI", "GBR", "CEU", "FIN",
                    "JPT", "CHS", "CHB", "CDX", "KHV",
                    "PJL", "STU", "BEB", "ITU", "GIH")

# Create integer labels for each SUPER-POPULATION
superpopulations_labels = c(1,1,1,1,1,1,1,
                            2,2,2,2,
                            3,3,3,3,3,
                            4,4,4,4,4,
                            5,5,5,5,5)
# Create integer labels for each SUB-POPULATION
sub_populations_labels = c(1,2,3,4,5,6,7,
                           21,22,23,24,
                           31,32,33,34,35,
                           41,42,43,44,45,
                           51,52,53,54,55)

popu_1kg = data.frame(cbind(superpopulations, 
                            sub_populations,
                            superpopulations_labels,
                            sub_populations_labels),
                      stringsAsFactors = F)

popu_1kg$sub_populations_labels = as.integer(popu_1kg$sub_populations_labels)
popu_1kg$superpopulations_labels = as.integer(popu_1kg$superpopulations_labels)

# Load list of platekeys within VCF fo run through fastPHASE
list_platekeys = read.table("./list_33714_platekeys.txt", stringsAsFactors = F)
list_platekeys = list_platekeys$V1
length(list_platekeys)
# 33714

# Load popu table
popu_table = read.csv("~/Documents/STRs/ANALYSIS/population_research/MAIN_ANCESTRY/GEL_60k_germline_dataset_fine_grained_population_assignment20200224.csv",
                      stringsAsFactors = F,
                      header = T)
dim(popu_table)
# 59464  36

# Dataframe with <PLATEKEY> and <BEST GUESS SUBPOPU>
df_platekey_popu = popu_table %>%
  filter(ID %in% list_platekeys) %>%
  select(ID, best_guess_predicted_ancstry)

# Merge with labels
df_platekey_popu = left_join(df_platekey_popu,
                             popu_1kg,
                             by = c("best_guess_predicted_ancstry" = "sub_populations"))
dim(df_platekey_popu)
# 33714  5

# Take the SUBPOPU labels for the list of platekeys KEEPING THE SAME ORDER!!!
index_platekeys = which(df_platekey_popu$ID == list_platekeys)

# Are identical?
identical(list_platekeys, df_platekey_popu$ID)
# FALSE 

identical(list_platekeys, df_platekey_popu$ID[index_platekeys])
# FALSE 

# Since it does not keep the same order I'll do it in cromagnon style
l_subpopus = c()
for(i in 1:length(list_platekeys)){
  l_subpopus = c(l_subpopus,
                 df_platekey_popu %>%
                   filter(ID %in% list_platekeys[i]) %>%
                   select(sub_populations_labels) %>%
                   pull())
}
length(l_subpopus)
# 33714

l_superpopus = c()
for(i in 1:length(list_platekeys)){
  l_superpopus = c(l_superpopus,
                 df_platekey_popu %>%
                   filter(ID %in% list_platekeys[i]) %>%
                   select(superpopulations_labels) %>%
                   pull())
}
length(l_superpopus)
# 33714

# Write into files
fileConn<-file("./60k_GRCH38_germline_mergedgVCF_chrX_67495316_67595385_unrelated_subpopus_labels.inp")
writeLines(as.character(l_subpopus), fileConn, sep = ' ')
close(fileConn)

fileConn<-file("./60k_GRCH38_germline_mergedgVCF_chrX_67495316_67595385_unrelated_superpopus_labels.inp")
writeLines(as.character(l_superpopus), fileConn, sep = ' ')
close(fileConn)

# Print the association of subpopulations and subpopulation labels
write.table(to_print, "./list_subpopulations_subpopulation_labels.txt", quote = F, row.names = F, col.names = F, sep = "\t")
