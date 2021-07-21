# Objective: new figure 4 that can replace old Table3
date ()
Sys.info ()[c("nodename", "user")]
commandArgs ()
rm (list = ls ())
R.version.string ## "R version 3.3.2 (2016-10-31)"

# libraries
library(dplyr); packageDescription ("dplyr", fields = "Version") #"0.8.5"
library(tidyverse); packageDescription ("tidyverse", fields = "Version") # "1.2.1
library(ggplot2); packageDescription ("ggplot2", fields = "Version") #"3.3.0"
library(reshape2)
library(tidyr)
library(scales)
library(cowplot)

# set the working directory
setwd("~/Documents/STRs/PAPERS/VALIDATION_PAPER/LANCET/APPEAL/")

# load number participant tested
df_tested = read.csv("./table_for_figure4_tested.tsv",
                     stringsAsFactors = F,
                     header = T,
                     sep = "\t")
dim(df_tested)
# 15  4

df_tested = df_tested %>% select(panel_id, disease, individuals_tested)
df_tested$disease = factor(df_tested$disease, levels = unique(df_tested$disease))

top = ggplot(df_tested, aes(x = disease, y = individuals_tested, fill = disease)) +
  geom_bar(stat = "identity", position=position_dodge()) +
  #geom_text(aes(label = individuals_tested), vjust = 1.5, colour = "black") +
  #scale_fill_manual(values=rep("darkgrey", length(df_tested$disease))) +
  scale_fill_manual(values=c("grey54","grey54","grey54","grey54","grey54","grey54","grey54","grey54","grey54","grey64","grey54","grey54","grey54","grey54","grey44")) +
  guides(fill = FALSE) +
  xlab("") + ylab("Individuals tested") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=0)) +
  scale_x_discrete(position = "top", labels = wrap_format(5)) 

# load number participant confirmed
df_confirmed = read.csv("./table_for_figure4_confirmed.tsv",
                     stringsAsFactors = F,
                     header = T,
                     sep = "\t")
dim(df_confirmed)
# 15  14

melt_confirmed = melt(df_confirmed)
melt_confirmed$disease = factor(melt_confirmed$disease, levels = unique(melt_confirmed$disease))

bottom = ggplot(melt_confirmed, aes(x = disease, label = value, y = value, fill = variable)) +
  geom_bar(stat = "identity") +
  guides(fill = FALSE) +
  xlab("") + ylab("") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=0)) +
  scale_y_reverse() +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) +
  ylab("Individuals confirmed")

png("new_Figure4.png",units="in", width=20, height=20, res=300)
print(plot_grid(top, bottom, ncol = 1))
dev.off()


