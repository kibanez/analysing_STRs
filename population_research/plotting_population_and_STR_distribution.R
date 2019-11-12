# Objective: take advantage from the work GEL bioinfo research group has done estimating ancestry for ~59,356 genomes
# and analyse and study the STR distribution across them
date ()
Sys.info ()[c("nodename", "user")]
commandArgs ()
rm (list = ls ())
R.version.string ## "R version 3.6.1 (2019-07-05)"

# libraries
library(dplyr); packageDescription ("dplyr", fields = "Version") #"0.8.3"
library(ggplot2); packageDescription ("ggplot2", fields = "Version") #"3.2.1"
library(reshape); packageDescription ("reshape", fields = "Version") #"0.8.8"
library(scatterplot3d); packageDescription("scatterplot3d", fields = "Version") # 0.3-41

# Set environment
setwd("/Users/kibanez/Documents/STRs/ANALYSIS/population_research/")

# Load population data
popu_table_enriched = read.csv("./population_info_enriched_59356_by_031019.tsv",
                               header = T,
                               sep = "\t",
                               stringsAsFactors = F)
dim(popu_table_enriched)
# 59356  20

# Load EH STR output data
# Research ~80K genomes, EH-v3.1.2- November 2019 -- but the ones that are ALSO included in the population info table
df = read.csv('~/Documents/STRs/ANALYSIS/population_research/EH_3.1.2_research_October2019_55419_genomes/merged/merged_population_genomes_53674_avg_EHv3.1.2.tsv',
              sep = '\t',
              stringsAsFactors = F,
              header = T)
dim(df)
# 2494 12

# Load thresholds
# STR annotation, threshold including the largest normal and the smallest pathogenic sizes reported
gene_annotation_normal = '/Users/kibanez/git/analysing_STRs/threshold_largest_normal_reported_research.txt'
gene_data_normal = read.table(gene_annotation_normal, stringsAsFactors=F, header = T)

gene_annotation_pathogenic = '/Users/kibanez/git/analysing_STRs/threshold_smallest_pathogenic_reported_research.txt'
gene_data_pathogenic = read.table(gene_annotation_pathogenic, stringsAsFactors=F, header = T)


# Functions
# Function that plots the STR repeat-size frequencies for a gene/locus across the cohort
plot_gene <- function(df_input, gene_name, gene_data_normal, gene_data_pathogenic, output_folder, assembly) {
  threshold_normal = gene_data_normal %>% filter(grepl(gene_name, locus)) %>% select(threshold) %>% unlist() %>% unname()
  threshold_pathogenic = gene_data_pathogenic %>% filter(grepl(gene_name, locus)) %>% select(threshold) %>% unlist() %>% unname()
  
  df_gene = df_input %>% filter(gene %in% gene_name)
  
  alt_number = df_gene$allele
  
  #df_gene_barplot = data.frame(number_repeats = alt_number, af = df_gene$AF)
  df_gene_barplot = data.frame(number_repeats = alt_number, af = df_gene$num_samples)
  
  # order by 'number of repetition'
  df_gene_barplot = unique(df_gene_barplot[order(df_gene_barplot[,1]),])
  
  rownames(df_gene_barplot) = df_gene_barplot$number_repeats
  
  df_gene_barplot$number_repeats = as.numeric(df_gene_barplot$number_repeats)
  
  gene_name = paste(gene_name, assembly, sep = '_')
  pdf_name = paste(output_folder, gene_name, sep = "/")
  pdf_name = paste(pdf_name, 'pdf', sep = ".")
  
  min_value = min(df_gene_barplot$number_repeats)
  max_value = max(threshold_pathogenic + 1, df_gene_barplot$number_repeats)
  
  aux_plot = ggplot(unique(df_gene_barplot), aes(x = number_repeats, y = af)) + 
    geom_bar(stat = "identity") + 
    geom_text(aes(label=af), vjust=0, size = 4, colour = "grey") +
    ylab("Allele frequency") + 
    xlab("Repeat sizes (repeat units)") + 
    ggtitle(gene_name) + 
    geom_vline(xintercept = threshold_normal, colour = 'blue', lty = 2) + 
    geom_vline(xintercept = threshold_pathogenic, colour = 'red', lty = 2) + 
    coord_cartesian(xlim = c(min_value,max_value))
  
  
  pdf(pdf_name)
  print(aux_plot)	
  dev.off()
  
}

# Let's first explore the data and see the population distribution across 59K
popu_table_enriched = popu_table_enriched %>%
  mutate(population = case_when(pred_african_ancestries >= 0.7 ~ "AFR",
                                pred_american_ancestries >= 0.7 ~ "AMR",
                                pred_european_ancestries >= 0.7 ~ "EUR",
                                pred_east_asian_ancestries >= 0.7 ~ "EAS",
                                pred_south_asian_ancestries >= 0.7 ~ "ASI",
                                pred_african_ancestries >= 0.25 & pred_american_ancestries >= 0.25  ~ "AFR-AMR",
                                pred_african_ancestries >= 0.25 & pred_european_ancestries >= 0.25 ~ "AFR-EUR",
                                pred_african_ancestries >= 0.25 & pred_east_asian_ancestries >= 0.25  ~ "AFR-EAS",
                                pred_african_ancestries >= 0.25 & pred_south_asian_ancestries >= 0.25  ~ "AFR-ASI",
                                pred_american_ancestries >= 0.25 & pred_european_ancestries >= 0.25  ~ "AMR-EUR",
                                pred_american_ancestries >= 0.25 & pred_east_asian_ancestries >= 0.25  ~ "AMR-EAS",
                                pred_american_ancestries >= 0.25 & pred_south_asian_ancestries >= 0.25  ~ "AMR-ASI",
                                pred_european_ancestries >= 0.25 & pred_american_ancestries >= 0.25  ~ "EUR-AMR",
                                pred_european_ancestries >= 0.25 & pred_east_asian_ancestries >= 0.25  ~ "EUR-EAS",
                                pred_european_ancestries >= 0.25 & pred_south_asian_ancestries >= 0.25  ~ "EUR-ASI",
                                pred_east_asian_ancestries >= 0.25 & pred_south_asian_ancestries >= 0.25 ~ "EAS-ASI"))

# With all population coctails
png("figures/population_distribution_56176_all_ancestries.png")
ggplot(data=popu_table_enriched %>% filter(!is.na(population)), 
       aes(x=pc2, y=pc1, colour = population)) +
  geom_hex(bins=300) +
  xlab("PC2 across 56,176 genomes") +
  ylab("PC1 across 56,176 genomes") +
  guides(fill = FALSE)
dev.off()

# Just "pure" ancestries
popu_table_enriched %>% filter(population %in% c("AFR", "EUR", "AMR", "EAS", "ASI")) %>% select(platekey) %>% unique() %>% pull() %>% length()
# 56176

png("figures/population_distribution_56176_pure_ancestries.png")
ggplot(data=popu_table_enriched %>% filter(population %in% c("AFR", "EUR", "AMR", "EAS", "ASI")) , 
       aes(x=pc2, y=pc1, colour = population)) +
  geom_hex(bins=300) +
  xlab("PC2 across 56,176 genomes") +
  ylab("PC1 across 56,176 genomes") +
  guides(fill = FALSE)
dev.off()

# Let's plot the raw numbers of each ancestry sub-cohort or sub-group
raw_numbers_popus = as.data.frame(table(popu_table_enriched$population))
colnames(raw_numbers_popus) = c("population", "Number of genomes")

raw_numbers_popus_pure = raw_numbers_popus %>% 
  filter(population %in% c("AFR", "EUR", "AMR", "EAS", "ASI")) 

# Pure ancestries
png("figures/barplot_pure_ancestry_groups_raw_numbers.png")
ggplot(raw_numbers_popus_pure, 
       aes(x = reorder(population, -`Number of genomes`), y = `Number of genomes`)) + 
  geom_bar(stat = "identity", aes(fill = population)) + 
  geom_text(aes(label=`Number of genomes`), vjust=-0.5, size = 4, colour = "grey") +
  ylab("Number of genomes") + 
  xlab("Ancestry cohorts - EHv3.1.2 - 56,176 total genomes") 
dev.off()

# All ancestries
png("figures/barplot_all_ancestry_groups_raw_numbers.png")
ggplot(raw_numbers_popus, 
       aes(x = reorder(population, -`Number of genomes`), y = `Number of genomes`)) + 
  geom_bar(stat = "identity", aes(fill = population)) + 
  geom_text(aes(label=`Number of genomes`), vjust=-0.5, size = 4, colour = "grey") +
  ylab("Number of genomes") + 
  xlab("Ancestry cohorts - EHv3.1.2 - 56,176 total genomes") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
dev.off()


# After plotting the ancestry distribution, let's see how different locus-STR-distribution change across different population cohorts
# let's take only the genomes for which we do have popu information 
l_genomes = popu_table_enriched %>%
  filter(population %in% c("AFR", "EUR", "AMR", "EAS", "ASI")) %>%
  select(platekey) %>%
  unique() %>%
  pull()

length(l_genomes)
# 55419

# Write into a file the list of platekeys for which we do have population info, in order to create the corresponding merge VCF files
write.table(l_genomes,
            file = "list_genomes_55419_pure_ancestry_info.txt",
            quote = F,
            row.names = F,
            col.names = F)

# Let's stratify the whole TSV in b37 and b38
df_b37 = df %>% filter(!grepl("chr", chr))
dim(df_b37)
# 1868 12

df_b38 = df %>% filter(grepl("chr", chr))
dim(df_b38)
# 2627 12

output_folder = "./figures/"

# ATN1 
plot_gene(df_b37, 'ATN1', gene_data_normal, gene_data_pathogenic, output_folder, "GRCh37")
plot_gene(df_b38, 'ATN1', gene_data_normal, gene_data_pathogenic, output_folder, "GRCh38")



