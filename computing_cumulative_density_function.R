date ()
Sys.info ()[c("nodename", "user")]
commandArgs ()
rm (list = ls ())
R.version.string ## "R version 3.3.2 (2016-10-31)"
library (dplyr); packageDescription ("dplyr", fields = "Version") #"0.5.0"
library(ggplot2); packageDescription ("ggplot2", fields = "Version") #"2.2.1"
library(grid); packageDescription ("grid", fields = "Version") #"3.3.2"
library(gridExtra); packageDescription ("gridExtra", fields = "Version") #"2.2.1"
library(reshape2); packageDescription ("reshape2", fields = "Version") #"1.4.2"
library(plyr); packageDescription ("plyr", fields = "Version") #"1.8.4"

options (width = 130)
#options (width = 240)

### FUNCTIONS


####

# empirical cumulative distribution function
# # D&C: analysis for 1 locus (FMR1)
# 
# fmr1_file = "/Users/kibanez/Documents/GEL_STR/GRCh38/17July2017/04082017_2200/tables/merged_5529genomes_july2017_not_normal_proband_consent_FMR1.tsv"
# fmr1_table = read.csv(fmr1_file, header = TRUE, sep = "\t")
# 
# # we need to retrieve info about the estimated expansion length (`alt_allele`) and compute the cumulative density
# 
# str_length_fmr1 = as.data.frame(fmr1_table %>% select(alt_allele) %>% pull())
# colnames(str_length_fmr1) = "str_fmr1"
# 
# # Melt the data frame
# ggdata <- melt(str_length_fmr1)
# 
# # Set the data frame, & add ecdf() data.
# ggdata <- ddply(ggdata, .(variable), transform, ecd=ecdf(value)(value))
# 
# # Create the CDF using ggplot.
# cdf <- ggplot(ggdata, aes(x=value)) + stat_ecdf(aes(colour=variable))
# 
# # Generate the CDF.
# cdf
# 
# ####

# empirical cumulative distribution function for all the loci

# Internal allele frequencies

# May2018
merged_tsv = '/Users/kibanez/Documents/GEL_STR/STR_internal_allele_frequencies/EH-offtarget-v2.5.5-duplicatesRemoved/merged/merged_loci_internal_STR_frequencies_may2018_20632_V4.vcf'

# August 2018
merged_tsv = '/Users/kibanez/Documents/GEL_STR/STR_internal_allele_frequencies/august2018/EH-offtarget-v2.5.5-duplicatesRemoved/merged/merged_loci_internal_STR_frequencies_august2018_30481_V4.tsv'

# January 2019
# Typo: I should have changed the month when merging the VCF files...
merged_tsv = '/Users/kibanez/Documents/GEL_STR/STR_internal_allele_frequencies/january2019/EH-offtarget-v2.5.5-duplicatesRemoved/merged/merged_loci_internal_STR_frequencies_may2018_42888_V4.tsv'

merged_table = read.csv(merged_tsv, header = TRUE, sep = '\t')

# Analyse the whole table by loci: split it into different genes
#path_output = "/Users/kibanez/Documents/GEL_STR/STR_internal_allele_frequencies/cumulativeFrequencies"
#path_output = "/Users/kibanez/Documents/GEL_STR/STR_internal_allele_frequencies/cumulativeFrequencies/august2018/"
path_output = "/Users/kibanez/Documents/GEL_STR/STR_internal_allele_frequencies/cumulativeFrequencies/january2019"

dir.create(path_output)

tables_output = paste(path_output, 'tables', sep = '/')

dir.create(tables_output)

normal_threshold_file = "/Users/kibanez/git/Rstudio/analysing_STRs/threshold_largest_normal_reported.txt"
pathogenic_threshold_file = "/Users/kibanez/git/Rstudio/analysing_STRs/threshold_smallest_pathogenic_reported.txt"

normal_threshold_table = read.csv(normal_threshold_file, header = TRUE, sep = "\t")
pathogenic_threshold_table = read.csv(pathogenic_threshold_file, header = TRUE, sep = "\t")

# List to all loci analysed
# We will only compute cumulative freqs for those we have defined normal and pathogenic thresholds
l_loci = merged_table %>% select(gene) %>% unique() %>% pull() %>% as.character()

for (i in 1:length(l_loci)){
  # freqs table
  locus_table = merged_table %>% filter(gene %in% l_loci[i])
  
  # locus name
  locus_name = as.character(locus_table$gene[1])

  # retrieve the threshold for each gene
  threshold_gene_normal = normal_threshold_table %>% filter(locus %in% locus_name) %>% select(threshold) %>% pull() %>% as.integer()
  threshold_gene_pathogenic = pathogenic_threshold_table %>% filter(locus %in% locus_name) %>% select(threshold) %>% pull() %>% as.integer()
  
  #colnames(locus_table) = c("chr","start","end","ref","alt","gene","repeat_number","repeat_motif","variant_type","ref_length_bp","repeat_alt_allele", "total_length_ref", "total_length_str","num_samples","af","list_samples")
  colnames(locus_table) = c("chr", "start", "end", "repeat_size", "gene", "ref", "alt", "Repeat_Motif", "num_samples", "AF", "list_samples")
  
  # Output name
  #pdf_output = paste(path_output, paste(locus_name, "_cumulativeFreqs_may2018_20632_V4.pdf", sep = ""), sep = "/")
  #pdf_output = paste(path_output, paste(locus_name, "_cumulativeFreqs_august2018_30481_V4.pdf", sep = ""), sep = "/")
  pdf_output = paste(path_output, paste(locus_name, "_cumulativeFreqs_january2019_42888_V4.pdf", sep = ""), sep = "/")
  
  # we need to retrieve info about the estimated expansion length (`alt_allele`) and compute the cumulative density
  str_length_locus = as.data.frame(locus_table %>% select(repeat_size) %>% pull())
  colnames(str_length_locus) = locus_name
  
  # Melt the data frame
  ggdata <- melt(str_length_locus)
  
  # Set the data frame, & add ecdf() data.
  ggdata <- ddply(ggdata, .(variable), transform, ecd=ecdf(value)(value))
  
  
  # Write down the cumulative freq absolute numbers
  
  table_file = paste(tables_output, paste(l_loci[i], '_cumFreqs.tsv' ,sep = ''), sep = '/')
  write.table(ggdata, file = table_file, row.names = F, col.names = T, sep = '\t', quote = F)
  
  # Create the CDF using ggplot.
  # kiki: compare what the value of `stat_ecdf` is en la esquinita superior
  cdf <- ggplot(ggdata, aes(x=value)) + stat_ecdf(aes(colour=variable)) + labs(x="Repeat length", y="Cumulative frequency") + geom_vline(xintercept = threshold_gene_normal, colour = 'blue', lty = 2) + geom_vline(xintercept = threshold_gene_pathogenic, colour = 'red', lty = 2)
  
  # Generate the CDF
  pdf(pdf_output)
  cdf
  dev.off()
  
}
