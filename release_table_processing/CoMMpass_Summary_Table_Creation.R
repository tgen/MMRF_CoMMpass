## This will build a full patient matrix table with primary result fields

# Load needed tools
library(tidyverse)
library(readxl)

seqFISH_tx_table <- "/Volumes/labs/MMRF/commpass/data_releases/MMRF_CoMMpass_IA21a/MMRF_CoMMpass_IA21a/summary_flat_files/MMRF_CoMMpass_IA21_genome_tumor_only_mm_igtx_pairoscope.tsv"
seqFISH_cna_table <- "/Volumes/labs/MMRF/commpass/data_releases/MMRF_CoMMpass_IA21a/MMRF_CoMMpass_IA21a/summary_flat_files/MMRF_CoMMpass_IA21_genome_gatk_cna_seqFISH.tsv"
qc_table <- "/Volumes/labs/MMRF/commpass/data_releases/MMRF_CoMMpass_IA21a/MMRF_CoMMpass_IA21a/README/MMRF_CoMMpass_IA21_Seq_QC_Summary.tsv"
curated_outcomes_table <- "/Volumes/labs/MMRF/commpass/data_releases/MMRF_CoMMpass_IA21a/MMRF_CoMMpass_IA21a/clinical_data_tables/CoMMpass_IA21_FlatFiles/MMRF_CoMMpass_IA21_Second_Line_Maint_Curated.tsv"

#### First build the full results table with all timepoints
#### Second subset a baseline table

#################################################
## Import translocation results select the columns of interest
ig_tx <- read_tsv(seqFISH_tx_table) %>% select(SAMPLE, ends_with("_CALL"), ends_with("_IGSOURCE"))

# Add column to count the number of detected translocations
ig_tx <- ig_tx %>% mutate(Ig_Translocation_Count = NSD2_CALL + CCND3_CALL + MYC_CALL + MAFA_CALL + CCND1_CALL + CCND2_CALL + MAF_CALL + MAFB_CALL)

# Split out Patient_ID and Visit_ID into seperate columns for merging events
ig_tx <- ig_tx %>% mutate(Patient_ID = str_sub(SAMPLE, start = 1, end = 9), 
                          Visit_ID = str_sub(SAMPLE, start = 1, end = 11))

# Add flag to indicate this sample has translocation data
ig_tx <- ig_tx %>% add_column("Has_Tx_Data" = 1)


#################################################
## Import copy number FISH table and select the columns of interest
cna_seqfish <- read_tsv(seqFISH_cna_table) %>% 
  select(SAMPLE, contains("Hyperdiploid"), contains("_1q21"), contains("13q14"), contains("13q34"), contains("17p13"), contains("CDKN2C"), contains("FAM46C"), contains("RB1"), contains("TP53"), contains("TRAF3")) %>% 
  select(-ends_with("_50percent")) %>% 
  rename_all( ~ str_replace(., "20percent", paste0("Call")))

# Split out Patient_ID and Visit_ID into seperate columns for merging events
cna_seqfish <- cna_seqfish %>% mutate(Patient_ID = str_sub(SAMPLE, start = 1, end = 9), 
                                      Visit_ID = str_sub(SAMPLE, start = 1, end = 11))

# Add flag to indicate this sample has WGS copy number data
cna_seqfish <- cna_seqfish %>% add_column("Has_CN_Data" = 1)


#################################################
## Import QC table to get reason for collection
qc_data <- read_tsv(qc_table, guess_max = 6000) %>% 
  select('Study Visit ID', Reason_For_Collection) %>% 
  rename("Visit_ID" = "Study Visit ID") %>% 
  unique()


#################################################
## Import curated outcomes data
outcomes <- read_tsv(curated_outcomes_table) %>% rename("Patient_ID" = "public_id") %>% add_column(Has_Outcome_Data = 1)

#################################################
## Merge tables into full table and add joint calculations when all data is available

# Merge seqFISH TX and CN tables
full_table <- full_join(ig_tx, cna_seqfish, by = c("Patient_ID", "Visit_ID", "SAMPLE"))

# Flag cytogenetic high risk samples
# Flag_0=Standard_risk_1=High_risk_(one_or_more_of_the_following_events:_del17p13_t(14;16)[MAF]_t(8;14)[MAFA]_t(14;20)[MAFB]_t(4;14)[WHSC1/MMSET/NSD2]
full_table <- full_table %>% mutate(Cytogenetic_High_Risk = case_when(Has_Tx_Data == 1 & Has_CN_Data == 1 & (SeqWGS_Cp_17p13_Call == 1 | NSD2_CALL == 1 | MAFA_CALL == 1 | MAF_CALL == 1 | MAFB_CALL == 1) ~ 1, 
                                                                      Has_Tx_Data == 1 & Has_CN_Data == 1 & SeqWGS_Cp_17p13_Call == 0 & NSD2_CALL == 0 & MAFA_CALL == 0 & MAF_CALL == 0 & MAFB_CALL == 0 ~ 0))


# Merge the outcomes data
full_table <- full_join(full_table, outcomes, by = "Patient_ID")

# Merge in reason for collection from QC table to filter to baseline visits
## WARNING: This should be a near final step to ensure all rows get a reason for collection added
full_table <- left_join(full_table, qc_data, by = "Visit_ID")


# Sort columns into useful order
full_table <- full_table %>% 
  select(Patient_ID, Visit_ID, SAMPLE, Reason_For_Collection, 
         starts_with("Has_"), 
         Ig_Translocation_Count, 
         Cytogenetic_High_Risk, 
         starts_with("NSD2"), 
         starts_with("CCND3"), 
         starts_with("MYC"), 
         starts_with("MAFA"), 
         starts_with("CCND1"), 
         starts_with("CCND2"), 
         starts_with("MAF"), 
         starts_with("MAFB"), 
         everything())


#################################################
## Create a table with just the baseline samples
baseline_table <- full_table %>% filter(Reason_For_Collection == "Baseline")


#################################################
## Write out tables to current working directory
write_tsv(full_table, "MMRF_CoMMpass_IA21_Summary_Table_All_Samples.tsv")
write_tsv(baseline_table, "MMRF_CoMMpass_IA21_Summary_Table_Baseline_Samples.tsv")
