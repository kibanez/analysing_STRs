# Objective: from EHv2.5.5 and EHv3.1.2 tables re STR calls, plot corresponding ROC curves

date ()
Sys.info ()[c("nodename", "user")]
commandArgs ()
rm (list = ls ())
R.version.string ## "R version 3.6.1 (2019-07-05)"

# libraries
library(dplyr)

# set the working directory
setwd("~/Documents/STRs/VALIDATION/raw_data/")

#The calculation has two steps:
# 1)Sort the observed outcomes by their predicted scores with the highest scores first
# 2)Calculate cumulative True Positive Rate (TPR) and True Negative Rate (TNR) for the ordered observed outcomes

# `labels` is a boolean vector with the actual classification of each case
# `scores` is a vector of real-valued prediction scores assigned by some classifier.
simple_roc <- function(labels, scores){
  labels <- labels[order(scores, decreasing=TRUE)]
  data.frame(TPR=cumsum(labels)/sum(labels), FPR=cumsum(!labels)/sum(!labels), labels)
}

