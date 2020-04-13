# Objective: select the best K, number of clusters, when running fastPHASE
date ()
Sys.info ()[c("nodename", "user")]
commandArgs ()
rm (list = ls ())
R.version.string ## "R version 3.6.1 (2019-07-05)"

require(dplyr); packageDescription ("dplyr", fields = "Version") #"0.8.3"
#On a local machine, the vignette can be accessed as follow:
browseVignettes("imputeqc")
system.file("extdata", package="imputeqc")

# Working directory
setwd("~/Documents/STRs/ANALYSIS/haplotyping/AR/")

# Workflow to select the best K


# 1- Generate a few test files with `make_test_files.R` which is enclosed to the package. 
# 5 test files each having 10% of masked genotypes are enough to start with.

# 2- Impute the missing genotypes in each test file with fastPHASE. Apply different K for each set of files.

# 3- Estimate the imputation quality with EstimateQuality() function and chose the K that minimizes the imputation error.
