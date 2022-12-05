#!/usr/bin/env Rscript --vanilla

# Load required modules
suppressPackageStartupMessages(require(tidyverse))
suppressPackageStartupMessages(require(optparse))

###############################################
## Define Options
###############################################

option_list = list(
  make_option(c("-m", "--mutation_counts"),
              type="character",
              default=NULL,
              help="File with count of mutations per gene per sample",
              metavar="filename"),
  make_option(c("-c", "--copy_number"),
              type="character",
              default=NULL,
              help="File with copy number per gene per sample",
              metavar="filename"),
  make_option(c("-e", "--gene_expression"),
              type="character",
              default=NULL,
              help="File with gene expression per gene per sample",
              metavar="filename"), 
  make_option(c("-d", "--gene_database"),
              type="character",
              default=NULL,
              help="File with gene expression per gene per sample",
              metavar="filename"), 
  make_option(c("-q", "--qc_table"),
              type="character",
              default=NULL,
              help="File molecular QC data (.xlsx)",
              metavar="filename"), 
  make_option(c("-s", "--subset"),
              type="character",
              default="Baseline",
              help="Visit Type to Plot (All, Baseline, Secondary) [Baseline]",
              metavar="Subset_of_Interest"), 
  make_option(c("-g", "--gene_name"),
              type="character",
              default=NULL,
              help="Official Gene Name to Query",
              metavar="Hugo_or_ENSG"), 
  make_option(c("-f", "--copyNumber_format"),
              type="character",
              default="Integer",
              help="How to plot copy number data (Log2 or Integer) [Integer] ",
              metavar="Log2_or_Integer"), 
  make_option(c("-r", "--rdata_object"),
              type="character",
              default="LOAD",
              help="Indicate if data should be saved to or loaded from existing Rdata",
              metavar="LOAD_or_SAVE"), 
  make_option(c("-n", "--rdata_name"),
              type="character",
              default=NULL,
              help="Prefix name of exising rData object or name to save data to",
              metavar="MyRdataFile")
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

###############################################
## Validate Inputs
###############################################


###############################################
## Read in input files or Load existing RData file
###############################################

if (opt$rdata_object == "SAVE") {
  # MMRF_CoMMpass_IA17_exome_vcfmerger2_All_Canonical_NS_Variants_ENSG_Mutation_Counts.tsv
  # Read in mutation file
  print("Reading Mutation Counts File...")
  mut <- read_tsv(opt$mutation_counts)
  
  # MMRF_CoMMpass_IA17_genome_gatk_cna_PerGene_LargestOverlap.tsv
  # Read in copy number file
  print("Reading Copy Number File...")
  cn <- read_tsv(opt$copy_number)
  
  # MMRF_CoMMpass_IA17_salmon_geneUnstranded_tpm.tsv
  # Read in expression file
  print("Reading Gene Expression File...")
  exp <- read_tsv(opt$gene_expression)
  
  # ensembl_v98_geneDB.tsv
  # Read in gene database file downloaded from biomart
  print("Reading Gene Database File...")
  gene_db <- read_tsv(opt$gene_database)
  
  # Read in QC table as this has the reason for collection
  print("Reading QC data File...")
  qc_data <- read_tsv(opt$qc_table, guess_max = 10000)
  qc_data <- qc_data %>% 
    rename("Patient_ID" = `KBase_Patient_ID`,
           "Study_Visit_ID" = `Study Visit ID`,
           "Reason_For_Collection" = `Reason_For_Collection`,
           "SampleName" = `QC Link SampleName`) %>%
    select(Patient_ID, Study_Visit_ID, Reason_For_Collection, SampleName) %>% 
    separate(SampleName, sep = "_", into = c("Study", "Patient", "Visit", "Source", "Fraction", "Drop_Tumor_Increment", "Drop_Tumor_Assay", "Drop_Library")) %>% 
    mutate(Subgroup = case_when(str_detect(Drop_Tumor_Increment, "T") ~ "Tumor", 
                                str_detect(Drop_Tumor_Increment, "C") ~ "Constitutional", 
                                TRUE ~ "Error")) %>% 
    filter(Subgroup == "Tumor") %>% 
    unite("Tumor_Specimen", Study:Fraction, sep = "_") %>% 
    select(-starts_with("Drop_")) %>% 
    unique()
  
  
  # Save the data files to an Rdata object, for fast recovery
  print("Saving Data as RData for fast access, in the future use (--rdata_object LOAD)...")
  save(mut, cn, exp, gene_db, qc_data, file = paste(opt$rdata_name, "RData", sep = "."))
} else {
  print("Loading Existing RData File...")
  load(paste(opt$rdata_name, "RData", sep = "."))
}


###############################################
## Set the gene to query
###############################################

#FGFR3
#gene_name <- "ENSG00000068078"
#CCND1
#gene_name <- "ENSG00000110092"
#TRAF3
#gene_name <- "ENSG00000131323"
#TP52
#gene_name <- "ENSG00000141510"
#MCL1
#gene_name <- "ENSG00000143384"

print(paste("Plotting Gene: ", opt$gene_name, sep = ""))

## Check if provided gene name is HGNC or ENSG
# Filter geneDB on HGNC symbol versus provided geneID
gene_check <- gene_db %>% filter(hgnc_symbol == opt$gene_name)
# Determine if the resulting table is empty or not
if (dim(gene_check)[1] == 0) {
  print("Provided Gene ID does not match HGNC Symbol in Database")
  print("Checking Against ENSG Identifiers")
  # Filter geneDB on ENSG ID versus provided geneID
  gene_check <- gene_db %>% filter(ensembl_gene_id == opt$gene_name)
  if (dim(gene_check)[1] == 0) {
    print("Provided Gene ID does not match Ensembl Gene ID in Database")
    print("ERROR - Exiting")
    exit()
  } else if (dim(gene_check)[1] == 1) {
    # Provided geneID was an ENSG ID, get associated HGNC symbol
    ensg_id <- opt$gene_name
    hgnc_id <- gene_check %>% pull(var = hgnc_symbol)
  } else {
    print("Provided Gene ID does not match Ensembl Gene ID or HGNC Gene ID in Database")
    print("ERROR - Exiting")
    q()
  }
} else if (dim(gene_check)[1] == 1) {
  # Provided geneID was a HGNC symbol, get associated ENSG ID
  ensg_id <- gene_check %>% pull(var = ensembl_gene_id)
  hgnc_id <- opt$gene_name
  print(paste("Provided HGNC Symbol Matches: ", ensg_id, sep = ""))
} else {
  print("Provided Gene ID does not match Ensembl Gene ID or HGNC Gene ID in Database")
  print("ERROR - Exiting")
  q()
}

plot_prefix <- paste(hgnc_id, ensg_id, opt$subset, sep = "_")

###############################################
## Process input files to select gene of interest
###############################################

## MUTATION COUNTS PROCESSING
# Select gene of interest from mutation table
gene_mut <- mut %>% filter(Gene == ensg_id)
print("Mut Filtered")

# Pivot the table, update gene name, process normal-tumor sampleID to get specimenID
gene_mut <- gene_mut %>% 
  pivot_longer(-Gene, names_to = "Tumor_Specimen", values_to = "Mutation_Count")

## COPY NUMBER PROCESSING
# Select gene of interest from the copy number table
gene_cn <- cn %>% filter(Gene == ensg_id)
print("CN Filtered")

# Pivot the table, process normal-tumor sampleID to get specimenID
gene_cn <- gene_cn %>% 
  pivot_longer(-Gene, names_to = "Tumor_Specimen", values_to = "Copy_Number_log2")

## GENE EXPRESSION PROCESSING
# Select gene of interest from the expression table
gene_exp <- exp %>% filter(Gene == ensg_id)
print("Exp Filtered")

# Pivot the table, process tumor sampleID to get specimenID
gene_exp <- gene_exp %>% 
  pivot_longer(-Gene, names_to = "Tumor_Specimen", values_to = "Expression_TPM")

###############################################
## Join the processed data tables into single table
###############################################

# Join data tables
gene_table <- full_join(gene_mut, gene_cn, by = "Tumor_Specimen")
gene_table <- full_join(gene_table, gene_exp, by = "Tumor_Specimen")
gene_table <- left_join(gene_table, qc_data, by = "Tumor_Specimen")

# Save data table
write_tsv(gene_table, paste(hgnc_id, ensg_id, "all_data.tsv", sep = "_"))

# Filter table
if (opt$subset == "All") {
  print("Plotting all visits with data")
} else if (opt$subset == "Baseline") {
  print("Plotting just the Baseline visits")
  gene_table <- gene_table %>% filter(Reason_For_Collection == "Baseline")
} else if (opt$subset == "Secondary") {
  print("Plotting just the secondary, non-baseline, visits")
  gene_table <- gene_table %>% filter(Reason_For_Collection != "Baseline")
}

# Sort joined table by mutation count
gene_table <- gene_table %>% arrange(Mutation_Count)

# Summarize table
gene_summary <- gene_table %>% 
  summarise(exp_mean = mean(Expression_TPM, na.rm=TRUE),
            exp_median = median(Expression_TPM, na.rm=TRUE),
            exp_min = min(Expression_TPM, na.rm=TRUE),
            exp_max = max(Expression_TPM, na.rm=TRUE),
            cn_mean = mean(Copy_Number_log2, na.rm=TRUE),
            cn_median = median(Copy_Number_log2, na.rm=TRUE),
            cn_min = min(Copy_Number_log2, na.rm=TRUE),
            cn_max = max(Copy_Number_log2, na.rm=TRUE))

# Capture the max expression value
maxEXPR <- gene_summary %>% pull(var = exp_max)
# Determine if it is above the minimum threshold
if ( maxEXPR < 128 ) {
  maxEXPR <- 128
}

# Capture the max copy number value
maxCN <- gene_summary %>% pull(var = cn_max)

# Capture the min copy number value
minCN <- gene_summary %>% pull(var = cn_min)

# Determin if the MAX CN is above the minimum thresholds
if ( maxCN < 1.321928095 ) {
  maxCNint <- 1.321928095
} else if ( maxCN >= 1.321928095) {
  maxCNint <- maxCN
}

if ( maxCN < 2.0 ) {
  maxCNlog2 <- 2.0
} else if ( maxCN >= 2.0) {
  maxCNlog2 <- maxCN
} 

# Determin if the MIN CN is below the minimum threshold
if ( minCN > -2.0 ) {
  minCN <- -2.0
}

###############################################
## Plot the Gene Expression Histogram
###############################################

ggplot(gene_exp, aes(log2(Expression_TPM+1))) + 
  geom_histogram(binwidth = 0.2) + 
  geom_vline(aes(xintercept=1.0), colour="#67a9cf", linetype="dashed", size = .5) + 
  annotate("text", x = 1.3, y = 30, label = "On", size = rel(2.5)) +
  annotate("text", x = 0.7, y = 30, label = "Off", size = rel(2.5)) +
  scale_x_continuous(limits=c(0, log2(maxEXPR + 1)), 
                     name = "Estimated Gene Expression Level [Log2(TPM+1)]") +
  scale_y_continuous(name = "Count") +
  ggtitle(label = plot_prefix) + 
  theme(plot.title = element_text(hjust = 0.5, size = 14), 
        axis.text = element_text(size=8, colour="black"), 
        axis.text.x = element_text(vjust= 0.5), 
        axis.title.y = element_text(size=10, vjust=1.5), 
        axis.title.x = element_text(size=10, vjust=-0.5), 
        legend.justification = "center", 
        legend.title = element_text(size=10), 
        legend.text = element_text(size=12),
        plot.margin=unit(c(.5,.5,.5,.5),"cm"))
ggsave(paste(plot_prefix, "expression_histogram.png", sep = "_"), 
       dpi = 150, 
       width = 100, 
       height = 100, 
       units = "mm")

###############################################
## Plot the Copy Number Histogram
###############################################
if ( opt$copyNumber_format == "Integer") {
  ggplot(gene_cn, aes((2 ^ Copy_Number_log2)*2)) + 
    geom_histogram(binwidth = 0.1, aes(y = ..density..)) + 
    geom_vline(aes(xintercept=1.788), colour="#ef8a62", linetype="dashed", size = .5) + 
    geom_vline(aes(xintercept=2.236), colour="#ef8a62", linetype="dashed", size = .5) + 
    #annotate("text", x = 0.2, y = 2.0, label = "Normal", size = rel(5.0)) + 
    annotate("text", x = 2.7, y = 1.75, label = "Gain", size = rel(2.5)) + 
    annotate("text", x = 1.3, y = 1.75, label = "Loss", size = rel(2.5)) + 
    scale_x_continuous(limits=c(0, ((2 ^ maxCNint)*2)), 
                       name = "Estimated Copies Per Cell") + 
    ggtitle(label = plot_prefix) + 
    theme(plot.title = element_text(hjust = 0.5, size = 14), 
          axis.text = element_text(size=8, colour="black"), 
          axis.text.x = element_text(vjust= 0.5), 
          axis.title.y = element_text(size=10, vjust=1.5), 
          axis.title.x = element_text(size=10, vjust=-0.5), 
          legend.justification = "center", 
          legend.title = element_text(size=10), 
          legend.text = element_text(size=12),
          plot.margin=unit(c(.5,.5,.5,.5),"cm"))
  
} else if ( opt$copyNumber_format == "Log2") {
  ggplot(gene_cn, aes(Copy_Number_log2)) + 
    geom_histogram(binwidth = 0.02)
}
ggsave(paste(plot_prefix, "copyNumber_histogram.png", sep = "_"), 
       dpi = 150, 
       width = 100, 
       height = 100, 
       units = "mm")

###############################################
## Plot the Three-way relationship
###############################################

# Remove rows with NA in the mutation column
gene_table <- gene_table %>% filter(Mutation_Count >= 0)

# Generate Custom Image depending on Copy Number Representation Requested
if ( opt$copyNumber_format == "Integer") {
  ggplot(gene_table, aes(log2(Expression_TPM + 1), (2 ^ Copy_Number_log2)*2)) +
    geom_point(aes(colour = factor(Mutation_Count)), size=3) +
    scale_colour_manual(values = c("grey", "#018571", "#8c510a","#d8b365", "#5ab4ac"), name="Mutation\nCount", limits=c(0,1,2,3,4)) +
    geom_hline(aes(yintercept=1.788), colour="#ef8a62", linetype="dashed", size = .5) +
    geom_hline(aes(yintercept=2.236), colour="#ef8a62", linetype="dashed", size = .5) +
    geom_vline(aes(xintercept=1.0), colour="#67a9cf", linetype="dashed", size = .5) + 
    annotate("text", x = 1.2, y = 0, label = "On", size = rel(5.0)) +
    annotate("text", x = 0.8, y = 0, label = "Off", size = rel(5.0)) +
    annotate("text", x = 0.2, y = 2.0, label = "Normal", size = rel(5.0)) +
    annotate("text", x = 0.2, y = 2.4, label = "Gain", size = rel(5.0)) +
    annotate("text", x = 0.2, y = 1.6, label = "Loss", size = rel(5.0)) + 
    scale_x_continuous(limits=c(0, log2(maxEXPR + 1)), 
                       name = "Estimated Gene Expression Level [Log2(TPM+1)]") + 
    scale_y_continuous(limits=c(0, ((2 ^ maxCNint)*2)), 
                       name = "Estimated Copies Per Cell") + 
    ggtitle(label = plot_prefix) + 
    theme(plot.title = element_text(hjust = 0.5, size = 20), 
          axis.text = element_text(size=14, colour="black"), 
          axis.text.x = element_text(vjust= 0.5), 
          axis.title.y = element_text(size=16, vjust=1.5), 
          axis.title.x = element_text(size=16, vjust=-0.5), 
          legend.justification = "center", 
          legend.title = element_text(size=10), 
          legend.text = element_text(size=12),
          plot.margin=unit(c(.5,0,.5,.5),"cm"))
} else if ( opt$copyNumber_format == "Log2") { 
  ggplot(gene_table, aes(log2(Expression_TPM + 1), Copy_Number_log2)) + 
    geom_point(aes(colour = factor(Mutation_Count)), size=3) + 
    scale_colour_manual(values = c("grey", "#018571", "#8c510a","#d8b365", "#5ab4ac"), name="Mutation\nCount", limits=c(0,1,2,3,4)) + 
    geom_hline(aes(yintercept=-0.1613), colour="#ef8a62", linetype="dashed", size = .5) + 
    geom_hline(aes(yintercept=0.1613), colour="#ef8a62", linetype="dashed", size = .5) + 
    geom_vline(aes(xintercept=1.0), colour="#67a9cf", linetype="dashed", size = .5) + 
    annotate("text", x = 1.2, y = -1.8, label = "On", size = rel(5.0)) + 
    annotate("text", x = 0.8, y = -1.8, label = "Off", size = rel(5.0)) + 
    annotate("text", x = 0.2, y = 0, label = "Normal", size = rel(5.0)) + 
    annotate("text", x = 0.2, y = 0.4, label = "Gain", size = rel(5.0)) + 
    annotate("text", x = 0.2, y = -0.4, label = "Loss", size = rel(5.0)) + 
    scale_x_continuous(limits=c(0, log2(maxEXPR + 1)), 
                       name = "Estimated Gene Expression Level [Log2(TPM+1)]") + 
    scale_y_continuous(limits=c(minCN, maxCNlog2), 
                       name = "Estimated Copy Number [Log2]") + 
    ggtitle(label = plot_prefix) + 
    theme(plot.title = element_text(hjust = 0.5, size = 20), 
          axis.text = element_text(size=14, colour="black"), 
          axis.text.x = element_text(vjust= 0.5), 
          axis.title.y = element_text(size=16, vjust=1.5), 
          axis.title.x = element_text(size=16, vjust=-.5), 
          legend.justification = "center", 
          legend.title = element_text(size=10), 
          legend.text = element_text(size=12))
} 

ggsave(paste(plot_prefix, "threeWay_plot.png", sep = "_"), 
       dpi = 300, 
       width = 324, 
       height = 184, 
       units = "mm")

