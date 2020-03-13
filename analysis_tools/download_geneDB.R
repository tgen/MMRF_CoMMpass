#!/usr/bin/env Rscript --vanilla

suppressPackageStartupMessages(require(biomaRt))

### To install biomart if not available
# if (!requireNamespace("BiocManager", quietly = TRUE))
#     install.packages("BiocManager")
# BiocManager::install("biomaRt")

## Download the gene database table from biomart
print("Downloading Gene Database File from BioMart...")
mart=useMart(biomart="ensembl",
             dataset="hsapiens_gene_ensembl",
             host = "http://sep2019.archive.ensembl.org")

gene_db <- getBM(attributes = c("ensembl_gene_id",
                           "ensembl_gene_id_version",
                           "chromosome_name",
                           "start_position",
                           "end_position",
                           "strand",
                           "band",
                           "external_gene_name",
                           "external_gene_source",
                           "gene_biotype",
                           "transcript_count",
                           "hgnc_id",
                           "hgnc_symbol",
                           "entrezgene_id"),
                 mart = mart)

write.table(gene_db, file = "ensembl_v98_geneDB.tsv", sep = "\t", row.names = FALSE)
print("Finished Writing Gene Database to Current Directory...")
