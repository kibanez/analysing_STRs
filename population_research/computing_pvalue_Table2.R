# Objective: run anova one way-like pvalue
# As we are working with percentages a beta regression sounds like a better strategy
date ()
Sys.info ()[c("nodename", "user")]
commandArgs ()
rm (list = ls ())
R.version.string ## "R version 3.6.3 (2020-02-29)"

# libraries
#install.packages("betareg")
#install.packages("emmeans")
#install.packages("GGally")
#install.packages("VGAM")

library(betareg); packageDescription ("betareg", fields = "Version") #3.1-4
library(emmeans); packageDescription ("emmeans", fields = "Version") #1.6.1
require(GGally)
require(VGAM)

# 
test = read.csv("~/Downloads/test_anova.tsv", stringsAsFactors = F, sep = "\t")
#Order factors by the order in data frame - otherwise, R will alphabetize them
test$superpopu = factor(test$superpopu, levels = unique(test$superpopu))


summary(m = vglm(percentage ~ superpopu, tobit(Upper = 800), data = test))





