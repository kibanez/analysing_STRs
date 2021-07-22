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

# Let's separate by panels
panel_a = df_tested %>% filter(panel_id %in% "A")
panel_b = df_tested %>% filter(panel_id %in% "B")
panel_c = df_tested %>% filter(panel_id %in% "C")
panel_d = df_tested %>% filter(panel_id %in% "D")


top = ggplot(df_tested, aes(x = disease, y = individuals_tested, fill = disease)) +
  geom_bar(stat = "identity", position=position_dodge()) +
  scale_fill_manual(values=c("grey54","grey54","grey54","grey54","grey54","grey54","grey54","grey54","grey54","grey64","grey54","grey54","grey54","grey54","grey44")) +
  guides(fill = FALSE) +
  xlab("") + ylab("Individuals tested") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=0)) +
  scale_x_discrete(position = "top", labels = wrap_format(5)) 

top_panel_a = ggplot(panel_a, aes(x = disease, y = individuals_tested, fill = disease)) +
  geom_bar(stat = "identity", position=position_dodge()) +
  scale_fill_manual(values=rep("grey54",length(panel_a$disease))) +
  guides(fill = FALSE) +
  xlab("") + ylab("") +
  #theme(axis.text.x = element_text(angle = 60, vjust = 0.5, hjust=.3)) +
  scale_x_discrete(labels = wrap_format(5))  +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) 

  
top_panel_b = ggplot(panel_b, aes(x = disease, y = individuals_tested, fill = disease)) +
  geom_bar(stat = "identity", position=position_dodge()) +
  scale_fill_manual(values=rep("grey54",length(panel_b$disease))) +
  guides(fill = FALSE) +
  xlab("") + ylab("") +
  theme(axis.text.x = element_text(angle = 60, vjust = 0.5, hjust=.3)) +
  scale_x_discrete(labels = wrap_format(5)) 

top_panel_c = ggplot(panel_c, aes(x = disease, y = individuals_tested, fill = disease)) +
  geom_bar(stat = "identity", position=position_dodge()) +
  scale_fill_manual(values=rep("grey54",length(panel_c$disease))) +
  guides(fill = FALSE) +
  xlab("") + ylab("") +
  theme(axis.text.x = element_text(angle = 60, vjust = 0.5, hjust=.3)) +
  scale_x_discrete(labels = wrap_format(5)) 

top_panel_d = ggplot(panel_d, aes(x = disease, y = individuals_tested, fill = disease)) +
  geom_bar(stat = "identity", position=position_dodge()) +
  scale_fill_manual(values=rep("grey54",length(panel_d$disease))) +
  guides(fill = FALSE) +
  xlab("") + ylab("") +
  theme(axis.text.x = element_text(angle = 60, vjust = 0.5, hjust=.3)) +
  scale_x_discrete(labels = wrap_format(5)) 

# load number participant confirmed
df_confirmed = read.csv("./table_for_figure4_confirmed.tsv",
                     stringsAsFactors = F,
                     header = T,
                     sep = "\t")
dim(df_confirmed)
# 15  14

panel_a_confirmed = df_confirmed %>% filter(disease %in% c("Hereditary ataxia", "Hereditary spastic paraplegia","Early onset and familial Parkinson's Disease",
                                                           "Complex Parkinsonism (includes pallido-pyramidal syndromes)",
                                                           "Early onset dystonia", "Early onset dementia",
                                                           "Amyotrophic lateral sclerosis or motor neuron disease", "Charcot-Marie-Tooth disease", "Ultra-rare undescribed monogenic disorders"))
panel_b_confirmed = df_confirmed %>% filter(grepl("Intellectual disability Plus", disease))
panel_c_confirmed = df_confirmed %>% filter(disease %in% c("Congenital myopathy", "Distal myopathies", "Congenital muscular dystrophy", "Skeletal muscle channelopathy"))
panel_d_confirmed = df_confirmed %>% filter(disease %in% "Intellectual disability")

melt_confirmed = melt(df_confirmed)
melt_confirmed$disease = factor(melt_confirmed$disease, levels = unique(melt_confirmed$disease))

melt_confirmed_a = melt(panel_a_confirmed)
melt_confirmed_b = melt(panel_b_confirmed)
melt_confirmed_c = melt(panel_c_confirmed)
melt_confirmed_d = melt(panel_d_confirmed)

melt_confirmed_a$disease = factor(melt_confirmed_a$disease, levels = unique(melt_confirmed_a$disease))
melt_confirmed_b$disease = factor(melt_confirmed_b$disease, levels = unique(melt_confirmed_b$disease))
melt_confirmed_c$disease = factor(melt_confirmed_c$disease, levels = unique(melt_confirmed_c$disease))
melt_confirmed_d$disease = factor(melt_confirmed_d$disease, levels = unique(melt_confirmed_d$disease))
  
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

plot_grid(top, bottom, ncol = 1)
ggsave("new_Figure4.svg")

# Split by panel ID
bottom_panel_a = ggplot(melt_confirmed_a, aes(x = disease, label = value, y = value, fill = variable)) +
  geom_bar(stat = "identity") +
  guides(fill = FALSE) +
  xlab("") + ylab("") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=0)) +
  scale_y_reverse() +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) 
  #ylab("Individuals confirmed")


bottom_panel_a
ggsave("new_Figure4_PanelA_bottom.svg")

top_panel_a
ggsave("new_Figure4_PanelA_top.svg")

png("new_Figure4_panelA_bottom.png",units="in", width=6, height=3, res=300)
print(bottom_panel_a)
dev.off()
png("new_Figure4_panelA_top.png",units="in", width=6, height=5, res=300)
print(top_panel_a)
dev.off()


bottom_panel_b = ggplot(melt_confirmed_b, aes(x = disease, label = value, y = value, fill = variable)) +
  geom_bar(stat = "identity") +
  guides(fill = FALSE) +
  xlab("") + ylab("") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=0)) +
  scale_y_reverse() +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) +
  ylab("Individuals confirmed")

bottom_panel_c = ggplot(melt_confirmed_c, aes(x = disease, label = value, y = value, fill = variable)) +
  geom_bar(stat = "identity") +
  guides(fill = FALSE) +
  xlab("") + ylab("") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=0)) +
  scale_y_reverse() +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) +
  ylab("Individuals confirmed")

bottom_panel_d = ggplot(melt_confirmed_d, aes(x = disease, label = value, y = value, fill = variable)) +
  geom_bar(stat = "identity") +
  guides(fill = FALSE) +
  xlab("") + ylab("") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=0)) +
  scale_y_reverse() +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) +
  ylab("Individuals confirmed")

# Merging top and bottom for each panel ID
print(plot_grid(top_panel_a, bottom_panel_a, ncol = 1))
print(plot_grid(top_panel_b, bottom_panel_b, ncol = 1))
print(plot_grid(top_panel_c, bottom_panel_c, ncol = 1))
print(plot_grid(top_panel_d, bottom_panel_d, ncol = 1))


png("new_Figure4_panelA.png",units="in", width=4, height=5, res=300)
print(plot_grid(top_panel_a, bottom_panel_a, ncol = 1))
dev.off()


plot_grid(top_panel_a, bottom_panel_a, ncol = 1)
ggsave("Figure4_panelA.svg")

plot_grid(top_panel_b, bottom_panel_b, ncol = 1)
ggsave("Figure4_panelB.svg")

plot_grid(top_panel_c, bottom_panel_c, ncol = 1)
ggsave("Figure4_panelC.svg")

plot_grid(top_panel_d, bottom_panel_d, ncol = 1)
ggsave("Figure4_panelD.svg")
