# Public Tools for Analysis of MMRF CoMMpass Data
These tools assume you have access to the MMRF Researcher Gateway (research.themmrf.org) 
and have downloaded the required files from the "Downloads" section.

### Three-Way Plot Generator
This tool will generate three-ways plots with gene expression on the x-axis, copy number on
the y-axis, and each data point is colored by the number of mutations observed in the gene.

Requirements: R with the following packages: tidyverse, optparse, biomaRt

#### Step 1 - Download Gene Database File from Ensembl BioMart
MMRF CoMMpass release IA17 and later use ensembl version 98 gene models, which are equivalent to 
gencode v32.
``` 
# From the command line in your analysis folder
RScript --vanilla download_geneDB.R
```
#### Step 2 - Read in Data Files Downloaded from the Researcher Gateway
``` 
Usage: RScript --vanilla three_way_plot.R [options]


Options:
	-m FILENAME, --mutation_counts=FILENAME
		File with count of mutations per gene per sample

	-c FILENAME, --copy_number=FILENAME
		File with copy number per gene per sample

	-e FILENAME, --gene_expression=FILENAME
		File with gene expression per gene per sample

	-d FILENAME, --gene_database=FILENAME
		File with gene expression per gene per sample

	-g HUGO_OR_ENSG, --gene_name=HUGO_OR_ENSG
		Official Gene Name to Query

	-f LOG2_OR_INTEGER, --copyNumber_format=LOG2_OR_INTEGER
		How to plot copy number data (Log2 or Integer) [integer] 

	-r LOAD_OR_SAVE, --rdata_object=LOAD_OR_SAVE
		Indicate if data should be saved to or loaded from existing Rdata

	-n MYRDATAFILE, --rdata_name=MYRDATAFILE
		Prefix name of exising rData object or name to save data to

	-h, --help
		Show this help message and exit
```
