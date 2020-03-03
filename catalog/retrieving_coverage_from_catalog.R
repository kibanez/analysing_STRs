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

con <- initOpencgaR(host="https://opencgainternal.gel.zone/opencga/", version="v1")
con <- opencgaLogin(opencga = con, userid = "mbleda", passwd = "",
                    autoRenew = TRUE, showToken = TRUE)
test <- fileClient(OpencgaR = con, action = "search", batch_size=5000,
                   params = list(study = "1000000032", name="LP3001783-DNA_A01.bam"))
test$stats$WHOLE_GENOME_COVERAGE$coverageSummary[[1]]