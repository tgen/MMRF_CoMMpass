# RNAseq Based Myeloma Purity Calculator
----

This is designed to determine if an RNAseq alignment file contains a pure monoclonal B cell population

### REQUIREMENTS

* samtools
* featurecounts (subread package)
* R
  * ggplot2
  * psych (for geometric mean calculation)
* Base Unix commands
  * awk
  * grep
  * sort

### Usage

```
  Usage: purityChecker_GRCh38_e98.sh [ options ]

  -h      Display Help

  Required for Sample Calculations  
  -t      Threads to use [1]  
  -b      Input RNAseq file from STAR alignment  
  -a      File type (BAM/CRAM) [BAM] - CRAM will be converted to BAM on the fly  
  -w      Workflow (Jetstream/Other) [Jetstream] - workflow type controls final location of outputs (./stats vs Current Directory)  
  -r      Resource Directory [defaults to executable path]  
  -d      Cleanup Temp files (Yes/No) [Yes] - Set to "No" to debug  

  Required to Build Needed Resource Files  
  -p      Prepare Resource Files (Yes/No) [No]  
  -g      Input GTF (only needed to prepare resource files)  
  -l      List of Genes not expressed in B-cells (only needed to prepare resource files)  
  -i      Immunoglobulin Loci in GTF format (only needed to prepare resource files)  
```  

#### Before the first run you must create a number of needed files
To do this you must have a GTF file (tested with ensembl v98 from JetStream Resources - Phoenix)
```
purityChecker_GRCh38_e98.sh -r /my/output/directory \
    -p Yes \
    -g geneModel.gtf \
    -l Non-Bcell_Contamination_GeneList_e98.txt \
    -i Immunoglobulin_GRCh38_Loci.gtf
```

#### For each RNAseq alignment file
``` 
purityChecker_GRCh38_e98.sh -t 10 \
    -b MySample_RNAseq_Alignment.bam \
    -a BAM \
    -w Other \
    -r /my/output/directory \
    -d Yes 
```
* CRAM files can be used as inputs but they are converted to a BAM on the fly as featureCounts cannot read a CRAM today
  * After generating the results the temp BAM file is deleted if the cleanup mode is set to "Yes"

### Concept

* Leverage the unique alignments in the immunoglobulin loci to separate monoclonal from polyclonal B-cell/plasma cell populations

We expect that a pure myeloma cell population should have a monoclonal expression of a single heavy and light chain molecule. We assume the unique heavy chain variable and constant gene segments will represent the vast majority of all RNAseq reads aligned to the IgH locus and similarly as single light chain variable element and constant region will represent the majority of the light chain alignments.  We expect that a heavy chain and light chain expressing patient will have single IgH VH and CH high frequency alleles and single high frequency VK or VL and CK or CL alleles as shown.

* Leverage the unique genes expressed in non B-cells that can be used to determine if the sample is contaminated with other unexpected white blood cells

We also calculate a "Non B-Cell" contamination index by calculating the RPKM (by hand using this tool) for a set of genes known to be present in other white blood cells.  We then calculate the geometric mean of these observed expression values to determine if the sample is likely contaminated with non-B cells.

* Provide an overview of individual samples to determine if they are pure monoclonal plasma cell populations devoid of non-B cell cellular contamination.

#### Expected Patterns from a High Purity anti-CD138 Sort from a Multiple Myeloma Patient
IgH locus showing most IGH variable aligned reads on a single VH segment and most IGH constant aligned reads on a single IgH constant gene
![IgH Locus](/myeloma_purityCalculator_RNAseq/images/Expected_Good_Myeloma_IgH.png)

IgL loci showing most IgK or IgL variable aligned reads on a single VK or VL segment and most IgK or IgL constant aligned reads on a single IgL constant gene
![IgL Loci](/myeloma_purityCalculator_RNAseq/images/Expected_Good_Myeloma_IgL.png)

#### Typical patterns seen in a low purity anti-CD138 sort from a multiple myeloma patient
IgH locus showing a poly-clonal representation with most IgH variable and constant regions having some signal, typically follows the normal B-cell expected distribuion
![IgH Locus](/myeloma_purityCalculator_RNAseq/images/Low_Purity_Myeloma_IgH.png)

IgL loci showing a poly-clonal represenation and the near-normal expectation of 66.6% kappa and 33.3% lambda
![IgL Loci](/myeloma_purityCalculator_RNAseq/images/Low_Purity_Myeloma_IgL.png)

### Notes

This script will process a bam file approximately every 5 minutes

To summarize results into a table use `purityChecker_Summary.sh` and to graph the results `purityChecker_Summary.R`

