
cc_100 = read.csv("~/Documents/STRs/ANALYSIS/population_research/PAPER/carriers/cc_pileup_100Kg/summary_cc_pileup_100Kg_30sept_VGD_KI.tsv", stringsAsFactors = F, header = T, sep = "\t")
dim(cc_100)

# Enrich with `is_unrel` data
l_unrel = read.table("~/Documents/STRs/ANALYSIS/population_research/MAIN_ANCESTRY/batch2/l_unrelated_55603_genomes_batch2.txt", stringsAsFactors = F)
l_unrel = l_unrel$V1
length(l_unrel)
# 55603

cc_100 = cc_100 %>% mutate(is_unrel = ifelse(platekey %in% l_unrel, "Yes", "No"))

# Enrich with popu data
popu_info = read.csv("~/Documents/STRs/ANALYSIS/population_research/MAIN_ANCESTRY/popu_merged_batch1_batch2_79849_genomes.tsv", 
                     stringsAsFactors = F, 
                     header = F, 
                     sep = "\t")

cc_100 = left_join(cc_100, popu_info, by = c("platekey" = "V1"))
colnames(cc_100)[15] = "popu"

batch2_genomes = popu_info %>% filter(V1 %in% l_unrel)

clin_data = read.csv("~/Documents/STRs/clinical_data/clinical_data/Main_RE_V10_and_Pilot_programmes.tsv", stringsAsFactors = F, header = T, sep = "\t")

# Disease_group info (and all other clinical characteristics) we've got for probands
# Let's take the familyIDs that have been recruited as Neuro in `disease_group`
l_fam_neuro = clin_data %>%
  filter(grepl("[Nn][Ee][Uu][Rr][Oo]", diseasegroup_list)) %>%
  select(rare_diseases_family_id) %>%
  unique() %>%
  pull()
length(l_fam_neuro)
# 14402

clin_data = clin_data %>% select(platekey, rare_diseases_family_id, diseasegroup_list)
clin_data = clin_data %>% 
  group_by(rare_diseases_family_id) %>% 
  mutate(is_neuro = ifelse(rare_diseases_family_id %in% l_fam_neuro, "Neuro", "NotNeuro")) %>% 
  ungroup() %>% 
  as.data.frame()
clin_data = unique(clin_data)
dim(clin_data)
# 109411  4

cc_100 = left_join(cc_100, 
                   clin_data %>% select(platekey, diseasegroup_list,is_neuro), 
                   by = "platekey")

write.table(cc_100, "~/Documents/STRs/ANALYSIS/population_research/PAPER/carriers/cc_pileup_100Kg/summary_cc_pileup_100Kg_30sept_VGD_KI_unrel.tsv", quote = F, col.names = T, row.names = F, sep = "\t")

# How many unrel genomes NOT NEURO do we have?
unrel_disease_group = clin_data %>% 
  filter(platekey %in% l_unrel, is_neuro %in% "NotNeuro")

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
# 29  3

unrel_htt = left_join(unrel_htt,
                      popu_info,
                      by = c("platekey" = "V1"))

# Write into a table
write.table(unrel_htt,
            "~/Documents/STRs/ANALYSIS/population_research/100K/carrier_freq/table_29_unrel_expanded_HTT_enriched_popu_diseasegroup.tsv",
            quote = F,
            row.names = F,
            col.names  = T,
            sep = "\t")
