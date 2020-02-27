# Objective: analyse repeat-sizes across other ancestries in the 1Kg genome project
date ()
Sys.info ()[c("nodename", "user")]
commandArgs ()
rm (list = ls ())
R.version.string ## "R version 3.6.1 (2019-07-05)"

# libraries
library(dplyr)

# Set working directory
setwd("~/Documents/STRs/ANALYSIS/population_research/1kg/")

# Load data
all_data = read.csv("1000G_2504_high_coverage.sequence.index.tsv",
                    sep = "\t",
                    header = T,
                    stringsAsFactors = F)
dim(all_data)
# 2504  22

table(all_data$POPULATION)
#ACB ASW BEB CDX CEU CHB CHS CLM ESN FIN GBR GIH GWD IBS ITU JPT KHV LWK MSL MXL PEL PJL PUR STU TSI YRI 
#96  61  86  93  99 103 105  94  99  99  91 103 113 107 102 104  99  99  85  64  85  96 104 102 107 108 

l_all_popus = all_data$POPULATION

# Let's enrich them with gender information
igsr_gender = read.csv("igsr_samples.tsv",
                       sep = "\t",
                       stringsAsFactors = F,
                       header = T)
dim(igsr_gender)
# 2504  9

all_data = left_join(all_data,
                     igsr_gender %>% select(Sample.name, Sex),
                     by = c("SAMPLE_NAME" = "Sample.name"))


for(i in 1:length(l_all_popus)){
  dataset = all_data %>%
    filter(POPULATION %in% l_all_popus[i])  %>%
    select(ENA_FILE_PATH, SAMPLE_NAME, POPULATION, Sex)
  
  number_samples = length(unique(dataset$SAMPLE_NAME))
  output_name = paste(paste(l_all_popus[i], as.character(number_samples), sep = "_1Kg_"), ".csv", sep = "samples")
  write.table(dataset, paste("./each_population/", output_name, sep = ""), quote = F, row.names = F, col.names = T, sep = ",")
}

