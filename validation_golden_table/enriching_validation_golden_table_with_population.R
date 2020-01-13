library("dplyr")
library("ggplot2")

val_data = read.csv("~/Documents/STRs/VALIDATION/STRVALIDATION_ALLDATA_2019-10-7_ALL_kibanez_enriched.tsv", stringsAsFactors = F, header = T, sep = "\t")
dim(val_data)
# 638  20

# Let's take the list of platekeys within the validation list
l_validation = unique(val_data$LP_Number)
length(l_validation)
# 256

popu_table = read.csv("~/Documents/STRs/ANALYSIS/population_research/population_and_super-population_definitions_across_59352_WGS_REv9_271119.tsv", stringsAsFactors = F, header = T, sep = "\t")
dim(popu_table)
# 59356  21

l_platekeys = val_data$LP_Number

l_platekeys[348] = "LP2000865-DNA_G07"

popu_table_subset = popu_table %>% select(platekey, population, pc1, pc2)

popu_table = left_join(val_data, popu_table_subset, by = c("LP_Number"="platekey"))
popu_table = popu_table %>% 
  filter(LP_Number %in% l_validation) %>% 
  select(LP_Number, population, pc1, pc2)
dim(popu_table)
# 638  4

# there are many of them that are NA, because they are Pilot
# Jan'2020 - we now have info for Pilot genomes
popu_pilot_table = read.csv("~/Documents/STRs/ANALYSIS/population_research/PILOT_ANCESTRY/FINE_GRAINED_RF_classifications_incl_superPOP_prediction_final20191216.csv",
                            stringsAsFactors = F,
                            sep = ",")

dim(popu_pilot_table)
# 4821  44

# Retrieve from here the df with platekey, bestGUESS_super_pop
popu_pilot_table = popu_pilot_table %>% 
  filter(ID %in% l_validation) %>%
  select(ID, bestGUESS_super_pop, PC_1, PC_2)
dim(popu_pilot_table)
# 63  4

colnames(popu_pilot_table) = c("LP_Number", "population", "pc1", "pc2")

# Before going beyond, keep only not NA rows in population in val_data_popu
index_notNA = which(!is.na(popu_table$population))
popu_table = popu_table[index_notNA,]
dim(popu_table)
# 154  4

val_data_popu = rbind(popu_table,
                      popu_pilot_table)
dim(val_data_popu)
# 217  4

val_data_popu = unique(val_data_popu)
dim(val_data_popu)
# 212  4

length(unique(val_data_popu$LP_Number))
# 212

# Unify or normalise all super populations nomenclature
table(val_data_popu$population)
#AFR     African         American AMR     AMR-EUR         ASI       EUR    European    South Asian 
#6          1           1      1          7               8          127          56        5

val_data_popu = val_data_popu %>% 
  mutate(merged_popu = case_when(population == "African" ~ "AFR",
                                 population == "European" ~ "EUR",
                                 population == "American" ~ "AMR",
                                 population == "South Asian" ~ "ASI",
                                 population == "East Asian" ~ "EAS",
                                 population == "AFR" ~ "AFR",
                                 population == "EUR" ~ "EUR",
                                 population == "AMR" ~ "AMR",
                                 population == "ASI" ~ "ASI",
                                 population == "AMR-EUR" ~ "AMR-EUR"))

write.table(val_data_popu, 
            "~/Documents/STRs/VALIDATION/population/STRVALIDATION_ALLDATA_2019-10-7_ALL_kibanez_enriched_with_popu_130120.tsv",
            sep = "\t", quote = F, row.names = F, col.names = T)


val_data_popu_distinct = distinct(val_data_popu, LP_Number, merged_popu)
dim(val_data_popu_distinct)
# 212  2

write.table(val_data_popu_distinct, 
            "~/Documents/STRs/VALIDATION/population/STRVALIDATION_ALLDATA_2019-10-7_ALL_kibanez_list_unique_platekeys_with_popu_130120.tsv",
            quote = F, col.names = T, row.names = F, sep = "\t")


# Let's plot them (for the ones we do have information)
# We need to plot them separately main and pilot, since PCs have been computed with diff number of PC (pilot with 6 PCs and main with 10 PCs)

val_data_popu_golden = val_data_popu %>% 
  filter(LP_Number %in% l_validation)


png("~/Documents/STRs/VALIDATION/population/population_distribution_154main_outof_217_all_ancestries.png")
ggplot(data=popu_table %>% filter(!is.na(population)), 
       aes(x=pc2, y=pc1, colour = population)) +
  #geom_hex(bins=100) +
  geom_point() +
  xlab("PC2 across 154 (main) out of 217 genomes") +
  ylab("PC1 across 154 (main) out of 217 genomes") +
  guides(fill = FALSE)
dev.off()

png("~/Documents/STRs/VALIDATION/population/population_distribution_63pilot_outof_217_all_ancestries.png")
ggplot(data=popu_pilot_table %>% filter(!is.na(population)), 
       aes(x=pc2, y=pc1, colour = population)) +
  #geom_hex(bins=100) +
  geom_point() +
  xlab("PC2 across 63 (pilot) out of 217 genomes") +
  ylab("PC1 across 63 (pilot) out of 217 genomes") +
  guides(fill = FALSE)
dev.off()

# What is happening for 44 genomes within the golden validation table? Whay we don't have info? they should be main....
l_genomes_with_popu = unique(c(popu_table$LP_Number, popu_pilot_table$LP_Number))
length(l_genomes_with_popu)
# 212

cuales = setdiff(unique(val_data$LP_Number),
                 l_genomes_with_popu)
length(cuales)
# 44

# Let's plot the raw numbers of each ancestry sub-cohort or sub-group
raw_numbers_popus = as.data.frame(table(val_data_popi$population))
colnames(raw_numbers_popus) = c("population", "Number of genomes")

# Pure ancestries
png("~/Documents/STRs/VALIDATION/population/barplot_population_groups_raw_numbers.png")
ggplot(raw_numbers_popus, 
       aes(x = reorder(population, -`Number of genomes`), y = `Number of genomes`)) + 
  geom_bar(stat = "identity", aes(fill = population)) + 
  geom_text(aes(label=`Number of genomes`), vjust=-0.5, size = 4, colour = "grey") +
  ylab("Number of genomes") + 
  xlab("149 out of 256 genomes with orthogonal validation") 
dev.off()

