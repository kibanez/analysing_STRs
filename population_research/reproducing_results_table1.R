
cc_100 = read.csv("~/Documents/STRs/ANALYSIS/population_research/PAPER/carriers/cc_pileup_100Kg/summary_cc_pileup_100Kg_30sept_VGD_KI.tsv", stringsAsFactors = F, header = T, sep = "\t")
dim(cc_100)

l_unrel = read.table("~/Documents/STRs/ANALYSIS/population_research/MAIN_ANCESTRY/batch2/l_unrelated_55603_genomes_batch2.txt", stringsAsFactors = F)

cc_100 = cc_100 %>% mutate(is_unrel = ifelse(platekey %in% l_unrel, "Yes", "No"))


#write.table(cc_100, "~/Documents/STRs/ANALYSIS/population_research/PAPER/carriers/cc_pileup_100Kg/summary_cc_pileup_100Kg_30sept_VGD_KI_unrel.tsv", quote = F, col.names = T, row.names = F, sep = "\t")


popu_info = read.csv("~/Documents/STRs/ANALYSIS/population_research/MAIN_ANCESTRY/popu_merged_batch1_batch2_79849_genomes.tsv", stringsAsFactors = F, header = F, sep = "\t")

cc_100 = left_join(cc_100, popu_info, by = c("platekey" = "V1"))

#write.table(cc_100, "~/Documents/STRs/ANALYSIS/population_research/PAPER/carriers/cc_pileup_100Kg/summary_cc_pileup_100Kg_30sept_VGD_KI_unrel.tsv", quote = F, col.names = T, row.names = F, sep = "\t")

batch2_genomes = popu_info %>% filter(V1 %in% l_unrel)

clin_data = read.csv("~/Documents/STRs/clinical_data/clinical_data/Main_RE_V10_and_Pilot_programmes.tsv", stringsAsFactors = F, header = T, sep = "\t")
clin_data = clin_data %>% select(platekey, diseasegroup_list)
clin_data = clin_data %>% group_by(platekey) %>% mutate(is_neuro = ifelse(grepl("[Nn][Ee][Uu][Rr][Oo]", diseasegroup_list), "Neuro", "NotNeuro")) %>% ungroup() %>% as.data.frame()
clin_data = unique(clin_data)
dim(clin_data)
# 109411  3

cc_100 = left_join(cc_100, clin_data, by = "platekey")
#write.table(cc_100, "~/Documents/STRs/ANALYSIS/population_research/PAPER/carriers/cc_pileup_100Kg/summary_cc_pileup_100Kg_30sept_VGD_KI_unrel.tsv", quote = F, col.names = T, row.names = F, sep = "\t")

unrel_disease_group = clin_data %>% filter(platekey %in% l_unrel)
unrel_disease_group = left_join(unrel_disease_group, batch2_genomes, by = c("platekey" = "V1"))
length(unique(unrel_disease_group$platekey))

write.table(unrel_disease_group, "~/Documents/STRs/ANALYSIS/population_research/PAPER/carriers/table_55603_unrel_genomes_enriched_popu_diseasegroup.tsv", sep = "\t", quote = F, row.names = F, col.names = F)

unrel_disease_group %>% filter(is_neuro %in% "NotNeuro", V2 %in% "AFR") %>% select(platekey) %>% unique() %>% pull() %>% length()
unrel_disease_group %>% filter(is_neuro %in% "NotNeuro", V2 %in% "AMR") %>% select(platekey) %>% unique() %>% pull() %>% length()
unrel_disease_group %>% filter(is_neuro %in% "NotNeuro", V2 %in% "EUR") %>% select(platekey) %>% unique() %>% pull() %>% length()
unrel_disease_group %>% filter(is_neuro %in% "NotNeuro", V2 %in% "EAS") %>% select(platekey) %>% unique() %>% pull() %>% length()
unrel_disease_group %>% filter(is_neuro %in% "NotNeuro", V2 %in% "SAS") %>% select(platekey) %>% unique() %>% pull() %>% length()


# Check unrel 29 genomes expanded in HTT
# Ask Ari in which table was done this
l_htt = read.table("~/Documents/STRs/ANALYSIS/population_research/100K/carrier_freq/list_29_expanded_after_QC_unrelated.tsv",
                   stringsAsFactors = F)
l_htt = l_htt$V1
length(l_htt)
# 29

# Enrich the list of HTT unrel genomes in 100kGP with popu and clinical info
unrel_htt = clin_data %>%
  filter(platekey %in% l_htt)
unrel_htt = unique(unrel_htt)
dim(unrel_htt)
# 1560  36