# The objective here is to enrich with clinical data all platekeys/participantIDs we do have in the PILOT cohort
date ()
Sys.info ()[c("nodename", "user")]
commandArgs ()
rm (list = ls ())
R.version.string ## "R version 3.3.2 (2016-10-31)"

# libraries
library(dplyr)

# setwd
setwd("~/Documents/STRs/clinical_data/pilot_clinical_data/")

# Load Pilot clinical data
data_platekeys = read.csv("~/Documents/STRs/clinical_data/pilot_clinical_data/data_freeze_Pilot_LK_RESEARCH/samples.csv",
                          sep = ",",
                          stringsAsFactors = FALSE,
                          header = TRUE)
dim(data_platekeys)
# 4833  2

reg_info = read.csv("~/Documents/STRs/clinical_data/pilot_clinical_data/data_freeze_Pilot_LK_RESEARCH/registration.csv",
                    sep = ",",
                    stringsAsFactors = FALSE,
                    header = TRUE)
dim(reg_info)
# 4877  23

reg_info = reg_info %>% select(gelId, sex, sex_at_birth, gender, relation_to_proband, isProband, gelFamilyId, disease_status, yearOfBirth)

diseases_info = read.csv("~/Documents/STRs/clinical_data/pilot_clinical_data/data_freeze_Pilot_LK_RESEARCH/diseases.csv",
                         sep = ",",
                         stringsAsFactors = FALSE,
                         header = TRUE)
dim(diseases_info)
# 2793  3

panel_info = read.csv("~/Documents/STRs/clinical_data/pilot_clinical_data/data_freeze_Pilot_LK_RESEARCH/panels.csv",
                      sep = ",",
                      stringsAsFactors = FALSE,
                      header = TRUE)
dim(panel_info)
# 3315  6

consent_info = read.csv("~/Documents/STRs/clinical_data/pilot_clinical_data/data_freeze_Pilot_LK_RESEARCH/pedigree.csv",
                      sep = ",",
                      stringsAsFactors = FALSE,
                      header = TRUE)
dim(consent_info)
# 17258  34

consent_info = consent_info %>% select(gelId, affectionStatus, review_outcome, qc_state, programmeConsent, consent_opinion)

## OUTPUT FILE
# First let's put all panel names as a list separated by ','
data_platekeys$gelID = as.character(data_platekeys$gelID)
to_write = left_join(data_platekeys,
                     reg_info %>% select(gelId, gelFamilyId, sex, sex_at_birth, gender, relation_to_proband, isProband, disease_status, yearOfBirth),
                     by = c("gelID" = "gelId"))
dim(to_write)
# 4833  10

# Let´s change all `relation_to_proband` == "" to "Proband"
to_write = to_write %>% group_by(gelID) %>% 
  mutate(biological_relation_to_proband = ifelse(isProband, "Proband", relation_to_proband)) %>% 
  as.data.frame()

to_write = to_write %>% select(gelID, gelFamilyId, plateKey, sex, biological_relation_to_proband, disease_status, yearOfBirth)
dim(to_write)
# 4833  7

diseases_info$gelID = as.character(diseases_info$gelID)
to_write = left_join(to_write,
                     diseases_info,
                     by = "gelID")
dim(to_write)
# 4974  9

consent_info$gelId = as.character(consent_info$gelId)
to_write = left_join(to_write,
                     consent_info %>% select(gelId, qc_state),
                     by = c("gelID" = "gelId"))
dim(to_write)
# 4974  10

panel_info$gelID = as.character(panel_info$gelID)
to_write = left_join(to_write,
                     panel_info %>% select(gelID, gelFamilyId, panel),
                     by = "gelID")
dim(to_write)
# 6202  12

# Put all panels separated by commas
list_panels = to_write %>% group_by(gelID) %>% summarise(panel_list = toString(panel)) %>% ungroup() %>% as.data.frame()
dim(list_panels)
# 4833  2

# Remove `panel_name` from to_write
to_write = to_write[,-12]
dim(to_write)
# 6202 11

to_write = left_join(to_write,
                     list_panels,
                     by = "gelID")

to_write = unique(to_write)

#plate_key.x	participant_id	programme	genome_build	programme_consent_status	rare_diseases_family_id	biological_relationship_to_proband	affection_status	participant_phenotypic_sex	year_of_birth	normalised_specific_disease	disease_sub_group	disease_group	panel_list
#LP2000950-DNA_H05	111000000	Rare Diseases	GRCh37	Consenting	111000000	Proband	Affected	Female	2005	Intellectual disability	Neurodevelopmental disorders	Neurology and neurodevelopmental disorders	Congenital hearing impairment (profound/severe), Congenital hearing impairment (profound/severe), Congenital hearing impairment (profound/severe), Congenital hearing impairment (profound/severe), Congenital hearing impairment (profound/severe), Congenital hearing impairment (profound/severe), Congenital hearing impairment (profound/severe), Congenital hearing impairment (profound/severe), Intellectual disability, Intellectual disability, Intellectual disability, Intellectual disability, Intellectual disability, Intellectual disability, Intellectual disability, Intellectual disability, Mitochondrial disorders, Mitochondrial disorders, Mitochondrial disorders, Mitochondrial disorders, Mitochondrial disorders, Mitochondrial disorders, Mitochondrial disorders, Mitochondrial disorders, Undiagnosed metabolic disorders, Undiagnosed metabolic disorders, Undiagnosed metabolic disorders, Undiagnosed metabolic disorders, Undiagnosed metabolic disorders, Undiagnosed metabolic disorders, Undiagnosed metabolic disorders, Undiagnosed metabolic disorders
to_write = to_write %>% select(gelID, plateKey, gelFamilyId.x, sex, biological_relation_to_proband, disease_status, yearOfBirth, specificDisease, ageOfOnset, qc_state, panel_list)

to_write %>% select(plateKey) %>% unique() %>% pull() %>% length()
# 4833

write.table(to_write, 
            file = "pilot_cohort_clinical_data_4833_genomes_withPanels_280919.tsv",
            quote = F,
            sep = "\t",
            row.names = F,
            col.names = T)

# Without panels
write.table(to_write[,-11], 
            file = "pilot_cohort_clinical_data_4833_genomes_removingPanels_280919.tsv",
            quote = F,
            sep = "\t",
            row.names = F,
            col.names = T)

