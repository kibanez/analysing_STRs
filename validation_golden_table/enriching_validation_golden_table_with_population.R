library("dplyr")
library("ggplot2")

val_data = read.csv("~/Documents/STRs/VALIDATION/STRVALIDATION_ALLDATA_2019-10-7_ALL_kibanez_enriched.tsv", stringsAsFactors = F, header = T, sep = "\t")
dim(val_data)

popu_table = read.csv("~/Documents/STRs/ANALYSIS/population_research/population_and_super-population_definitions_across_59352_WGS_REv9_271119.tsv", stringsAsFactors = F, header = T, sep = "\t")
dim(popu_table)

l_platekeys = val_data$LP_Number

l_platekeys[348] = "LP2000865-DNA_G07"

popu_table_subset = popu_table %>% select(platekey, population, pc1, pc2)

val_data_popi = left_join(val_data, popu_table_subset, by = "platekey")

write.table(val_data_popi, "~/Documents/STRs/VALIDATION/population/STRVALIDATION_ALLDATA_2019-10-7_ALL_kibanez_enriched_with_popu.tsv", sep = "\t", quote = F, row.names = F, col.names = T)


val_data_popi_distinct = distinct(val_data_popi, platekey, population)

write.table(val_data_popi_distinct, "~/Documents/STRs/VALIDATION/population/STRVALIDATION_ALLDATA_2019-10-7_ALL_kibanez_list_unique_platekeys_with_popu.tsv", quote = F, col.names = T, row.names = F, sep = "\t")


# Let's plot them (for the ones we do have information)
png("~/Documents/STRs/VALIDATION/population/population_distribution_149_outof_256_all_ancestries.png")
ggplot(data=val_data_popi %>% filter(!is.na(population)), 
       aes(x=pc2, y=pc1, colour = population)) +
  #geom_hex(bins=100) +
  geom_point() +
  xlab("PC2 across 149 out of 256 genomes") +
  ylab("PC1 across 149 out of 256 genomes") +
  guides(fill = FALSE)
dev.off()


# Let's plot the raw numbers of each ancestry sub-cohort or sub-group
raw_numbers_popus = as.data.frame(table(val_data_popi$population))
colnames(raw_numbers_popus) = c("population", "Number of genomes")

# Pure ancestries
png("barplot_population_groups_raw_numbers.png")
ggplot(raw_numbers_popus, 
       aes(x = reorder(population, -`Number of genomes`), y = `Number of genomes`)) + 
  geom_bar(stat = "identity", aes(fill = population)) + 
  geom_text(aes(label=`Number of genomes`), vjust=-0.5, size = 4, colour = "grey") +
  ylab("Number of genomes") + 
  xlab("149 out of 256 genomes with orthogonal validation") 
dev.off()

