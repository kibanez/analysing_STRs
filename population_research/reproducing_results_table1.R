# Objective: compute total number of unrel genomes, unrel not neuro genomes
# also, by ethnicity
date ()
Sys.info ()[c("nodename", "user")]
commandArgs ()
rm (list = ls ())
R.version.string ## "R version 3.6.1 (2019-07-05)"

cc_100 = read.csv("~/Documents/STRs/ANALYSIS/population_research/PAPER/carriers/cc_pileup_100Kg/summary_cc_pileup_100Kg_30sept_VGD_KI.tsv", stringsAsFactors = F, header = T, sep = "\t")
dim(cc_100)
# 1783  13

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
clin_data = unique(clin_data)
dim(clin_data)
# 109411  3

# Load platekey-pid-famID table we created to fish platekeys not included in further RE releases
clin_metadata = read.csv("~/Documents/STRs/clinical_data/clinical_data/merged_RE_releases_and_Pilot_PID_FID_platekey.tsv",
                         stringsAsFactors = F,
                         sep = "\t",
                         header = T)
dim(clin_metadata)
# 621704  4

# Include or enrich `clin_data` with extra platekeys, to associate platekey <-> famID
clin_data = full_join(clin_data,
                      clin_metadata %>% select(platekey, participant_id, rare_diseases_family_id),
                      by = "platekey")
clin_data = unique(clin_data)
dim(clin_data)
# 149776  5

# First let's unite `rare_diseases_family_id` columns into 1
clin_data = clin_data %>%
  group_by(rare_diseases_family_id.x) %>%
  mutate(famID = ifelse(is.na(rare_diseases_family_id.x), rare_diseases_family_id.y, rare_diseases_family_id.x)) %>%
  ungroup() %>%
  as.data.frame()

# Now we've got complete famID, let's define whether each platkey is neuro or not
clin_data = clin_data %>% 
  group_by(famID) %>% 
  mutate(is_neuro = ifelse(famID %in% l_fam_neuro, "Neuro", "NotNeuro")) %>% 
  ungroup() %>% 
  as.data.frame()

cc_100 = left_join(cc_100, 
                   clin_data %>% select(platekey, diseasegroup_list,is_neuro), 
                   by = "platekey")

write.table(cc_100, "~/Documents/STRs/ANALYSIS/population_research/PAPER/carriers/cc_pileup_100Kg/summary_cc_pileup_100Kg_30sept_VGD_KI_unrel.tsv", quote = F, col.names = T, row.names = F, sep = "\t")

# How many unrel genomes NOT NEURO do we have?
clin_data %>% filter(platekey %in% l_unrel, is_neuro %in% "NotNeuro") %>% select(platekey) %>% unique() %>% pull() %>% length()
# 37888
clin_data %>% filter(platekey %in% l_unrel, is_neuro %in% "Neuro") %>% select(platekey) %>% unique() %>% pull() %>% length()
# 17715

# Create an unrel (with 55603 genomes) clin_data table
unrel_disease_group = clin_data %>%
  filter(platekey %in% l_unrel) %>%
  select(platekey, famID, diseasegroup_list, is_neuro)
unrel_disease_group = unique(unrel_disease_group)
dim(unrel_disease_group)
# 55603  4

# Enrich with popu
unrel_disease_group = left_join(unrel_disease_group, batch2_genomes, by = c("platekey" = "V1"))
colnames(unrel_disease_group)[5] = "popu"
length(unique(unrel_disease_group$platekey))
# 55603

write.table(unrel_disease_group, "~/Documents/STRs/ANALYSIS/population_research/PAPER/carriers/table_55603_unrel_genomes_enriched_popu_diseasegroup.tsv", sep = "\t", quote = F, row.names = F, col.names = F)

unrel_disease_group %>% filter(is_neuro %in% "NotNeuro", popu %in% "AFR") %>% select(platekey) %>% unique() %>% pull() %>% length()
# 1295
unrel_disease_group %>% filter(is_neuro %in% "NotNeuro", popu %in% "AMR") %>% select(platekey) %>% unique() %>% pull() %>% length()
# 704
unrel_disease_group %>% filter(is_neuro %in% "NotNeuro", popu %in% "EUR") %>% select(platekey) %>% unique() %>% pull() %>% length()
# 32098
unrel_disease_group %>% filter(is_neuro %in% "NotNeuro", popu %in% "EAS") %>% select(platekey) %>% unique() %>% pull() %>% length()
# 314
unrel_disease_group %>% filter(is_neuro %in% "NotNeuro", popu %in% "SAS") %>% select(platekey) %>% unique() %>% pull() %>% length()
# 2947

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
