# This script will manipulate and summarize the CoMMpass Release File Manifest available in IA17+ releases

# Load required R libraries
library(tidyverse)

# Manually define the name of the manifest file
manifest_file <- "MMRF_CoMMpass_IA19_package_file_manifest.tsv"

# Read the manifest file into memory
data <- read_tsv(manifest_file)

# Make simple file for perVisit summaries (remove duplicate lines created from multiple data formats extracted from individual assays)
simple <- data %>% 
  select(project, visit, rgsm, subgroup, gltype) %>% 
  separate(rgsm, into = c("Study", "Patient", "Visit", "SampleSource", "Fraction", "Increment"), sep = "_") %>% 
  select(-Study, -Patient, -Visit) %>% 
  unique()

# get counts for unique patients and visits
patient_count <- simple %>% select(project) %>% unique() %>% count()
visit_count <- simple %>% select(visit) %>% unique() %>% count()


# Write out table for excel pivot summary
write_tsv(simple, "Release_Summary.tsv")
