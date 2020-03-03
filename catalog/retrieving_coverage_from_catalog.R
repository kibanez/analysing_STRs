# Objective: retrieve coverage info for platekeys from catalog
date ()
Sys.info ()[c("nodename", "user")]
commandArgs ()
rm (list = ls ())
R.version.string ## "R version 3.6.1 (2019-07-05)"

# installed the following packages
#install.packages("configr")
#install.packages("~/Downloads/opencgaR_1.4.0.tar.gz", repos = NULL, type = 'source')

# libraries
library(opencgaR)
library(dplyr)
library(tidyr)
library(purrr)

# working directory
setwd("~/Documents/STRs/VALIDATION/")

con <- initOpencgaR(host="https://opencgainternal.gel.zone/opencga/", version="v1")
con <- opencgaLogin(opencga = con, userid = "interpretationpipeline", passwd = "kx6SHa26",
                    autoRenew = TRUE, showToken = TRUE)

# load GEL validation platekeys
l_val_data = read.csv("./list_all_unique_254_genomes_GEL.txt", stringsAsFactors = F, header = F)
l_val_data = unique(l_val_data$V1)
length(l_val_data)
# 254

#write.table(l_val_data, "./list_all_unique_255_genomes_GEL.txt", quote = F, row.names = F, col.names = F)

df_coverage = c()
for (i in 1:length(l_val_data)){
  platekey = paste(l_val_data[i], ".bam", sep = "")
  
  # study '32 == RD GRCh38 dataset
  test <- fileClient(OpencgaR = con, action = "search", batch_size=5000,
                     params = list(study = "1000000032", name=platekey))
  if (!is.null(test)){
    df_coverage = test$stats$WHOLE_GENOME_COVERAGE$coverageSummary[[1]]
    autosome_cov = df_coverage %>% filter(scope %in% "autosomes") %>% select(avg) %>% pull()
    df_coverage = rbind(df_coverage,
                        cbind(platekey, autosome_cov))
    
  }
}
