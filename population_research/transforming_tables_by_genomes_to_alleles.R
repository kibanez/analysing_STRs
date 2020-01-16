# Objective: each row in the merged tables is a `repeat-size` together with the list of genomes having that repeat-size in an allele
# We do want to transform the table: First, splitting by locus, and then by having an allele per row
# Therefore, analysing an autosomal gene, we would expect to have 2 times alleles (and rows) for X genomes

# Also, important to check for duplicates. We shouldn't have, the input data is coming from the population aggregated file Loukas' group has been working on.

date ()
Sys.info ()[c("nodename", "user")]
commandArgs ()
rm (list = ls ())
R.version.string ## "R version 3.6.1 (2019-07-05)"

# libraries
library(dplyr)

# Functions
source("~/git/analysing_STRs/functions/transforming_genomeRow_to_alleleRow.R")

# Set working directory
#working_dir="~/Documents/STRs/ANALYSIS/population_research/EH_3.1.2_research_October2019_55419_genomes/unrelated_probands_and_cancer/AFR/"
#working_dir="~/Documents/STRs/ANALYSIS/population_research/EH_3.1.2_research_October2019_55419_genomes/unrelated_probands_and_cancer/AMR/"
#working_dir="~/Documents/STRs/ANALYSIS/population_research/EH_3.1.2_research_October2019_55419_genomes/unrelated_probands_and_cancer/EAS/"
#working_dir="~/Documents/STRs/ANALYSIS/population_research/EH_3.1.2_research_October2019_55419_genomes/unrelated_probands_and_cancer/ASI/"
working_dir="~/Documents/STRs/ANALYSIS/population_research/EH_3.1.2_research_October2019_55419_genomes/unrelated_probands_and_cancer/EUR/"
setwd(working_dir)

# Load data
#merged_file = "~/Documents/STRs/ANALYSIS/population_research/EH_3.1.2_research_October2019_55419_genomes/unrelated_probands_and_cancer/AFR/merged/merged_population_genomes_unrelated_probands_and_cancer_1136_avg_EHv3.1.2_AFR.tsv"
#merged_file = "~/Documents/STRs/ANALYSIS/population_research/EH_3.1.2_research_October2019_55419_genomes/unrelated_probands_and_cancer/AMR/merged/merged_population_genomes_unrelated_probands_and_cancer_384_avg_EHv3.1.2_AMR.tsv"
#merged_file = "~/Documents/STRs/ANALYSIS/population_research/EH_3.1.2_research_October2019_55419_genomes/unrelated_probands_and_cancer/EAS/merged/merged_population_genomes_unrelated_probands_and_cancer_227_avg_EHv3.1.2_EAS.tsv"
#merged_file = "~/Documents/STRs/ANALYSIS/population_research/EH_3.1.2_research_October2019_55419_genomes/unrelated_probands_and_cancer/ASI/merged/merged_population_genomes_unrelated_probands_and_cancer_2678_avg_EHv3.1.2_ASI.tsv"
merged_file = "~/Documents/STRs/ANALYSIS/population_research/EH_3.1.2_research_October2019_55419_genomes/unrelated_probands_and_cancer/EUR/merged/merged_population_genomes_unrelated_probands_and_cancer_26033_avg_EHv3.1.2_EUR.tsv"

# Transform genome-row tables to allele-row tables
transform_genome_to_allele_table(working_dir,merged_file)

list_files_tables = sort(list.files(pattern = "\\_simplified.tsv$"))

# Once the above function creates for each locus a table having an allele info per row, we will proceed to check for duplicates 
for (i in 1:length(list_files_tables)){
  merged_data = read.csv(list_files_tables[i],
                         sep = '\t',
                         header = T,
                         stringsAsFactors = F)


  # Keep with Cancer and RD programmes, filter out `.`
  merged_data = merged_data %>% filter(!programme %in% ".")
 
  output_file = paste(strsplit(list_files_tables[i], ".tsv")[[1]], "cancer_and_RD.tsv", sep = "_")
  write.table(merged_data,
              output_file,
              quote = F,
              row.names = F,
              col.names = T,
              sep = "\t") 
}
