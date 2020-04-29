# Objective: plot a bubble plot or 1:1 plot to analyse Hiseq vs Novaseq performance from STRs point of view
date ()
Sys.info ()[c("nodename", "user")]
commandArgs ()
rm (list = ls ())
R.version.string ## "R version 3.6.3(2020-02-29)"

library(ggplot2); packageDescription ("ggplot2", fields = "Version") #"2.2.1"
library(grid); packageDescription ("grid", fields = "Version") #"3.3.2"
library(gridExtra); packageDescription ("gridExtra", fields = "Version") #"2.3"
library(reshape2); packageDescription ("reshape2", fields = "Version") #"1.4.3"
require(dplyr); packageDescription ("dplyr", fields = "Version") #"0.7.4"
library(lattice); packageDescription ("lattice", fields = "Version") #"0.20-35"
library(tidyverse); packageDescription ("tidyverse", fields = "Version") #"1.2.1"
library(RColorBrewer); packageDescription ("RColorBrewer", fields = "Version") #"1.1-2"
library(cowplot); packageDescription ("cowplot", fields = "Version") #"1.0.0"

# Set working environment
setwd("/Users/kibanez/Documents/STRs/ANALYSIS/novaseq/STRs/")

# Load data

# Hiseq data
ash_hiseq = read.csv("./HG002_Edico.HiSeq.repeats.tsv",
                     stringsAsFactors = F, 
                     sep = "\t",
                     header = F)

pg_hiseq = read.csv("./NA12878_S1_Edico.HiSeq.repeats.tsv",
                    stringsAsFactors = F,
                    sep = "\t",
                    header = F)

# Novaseq1
ash_nova1 = read.csv("./LP4100018-DNA_E11_Edico.Ash.NovaSeq.repeats.tsv",
                     stringsAsFactors = F,
                     sep = "\t",
                     header = F)

pg_nova1= read.csv("./LP4100018-DNA_E05_Edico.NovaSeq.repeats.tsv",
                    stringsAsFactors = F, 
                    sep = "\t",
                    header = F)



# Novaseq2
ash_nova2 = read.csv("./LP4100018-DNA_F11_Edico.Ash.NovaSeq.repeats.tsv",
                     stringsAsFactors = F,
                     sep = "\t",
                     header = F)

pg_nova2= read.csv("./LP4100016-DNA_D02_Edico.NovaSeq.repeats.tsv",
                    stringsAsFactors = F, 
                    sep = "\t",
                    header = F)


merged_all = rbind(ash_hiseq,
                   pg_hiseq,
                   ash_nova1,
                   pg_nova1,
                   ash_nova2,
                   pg_nova2)
