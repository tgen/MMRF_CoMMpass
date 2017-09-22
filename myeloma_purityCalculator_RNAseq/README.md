# RNAseq Based Myeloma Purity Calculator
----

This is designed to determine if an RNAseq result bam file contains a pure monoclonal B cell population

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

* Ensure all variables at the top of the script have the path variables set correctly
  * RESULTS_TABLE=/YourPath/purityChecker_tabulatedResults.txt
  * BUILD_DIR=/YourPath/MMRF_CoMMpass/myelomaPurityChecker/
    * This is the folder where we expect a number of files to limit the number of required pointer variables
* On the first run you must create a number of needed files
  * To do this you must have a GTF file, use the related package to create the MMRF specific reference genome and GTF
  * Then set the build variable to create the needed files and ensure you point to the correct GTF
    * BUILD_FILES=Yes
    * GTF=/YourPath/Homo_sapiens.GRCh37.74.gtf.hs37d5.EGFRvIII.gtf
* To test a sample the script must be called in a directory with a RNAseq BAM file (This should be aligner agnastic)
  * Make sure the build variable is unset to ensure the required files are not created again
    * BUILD_FILES=No
  * Make sure the loop function search parameter is looking for the proper file extension matching your BAM files
    * for INPUT_BAM in `ls *.starAligned.final.bam`
    * REPLACE .starAligned.final.bam with the string that matches the ending of your BAM files
* The script will loop over in serial all bam files in the folder.  At TGen we submit each bam as independent jobs on our cluster but this is designed as a stand alone that should work in any enviroment.

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
![IgL Loci](/myeloma_purityCalculator_RNAseq/images/Low_Purity_Good_Myeloma_IgL.png)

### Notes

All other required files for GRCh37 are provided but you can build the required annotation files by changing the BUILD_FILES flag to "YES"

Currently, the script will process all .bam files in a folder from which the script is called

This script will process a bam file approximately every 5 minutes

If you need to graph the summary table you can use the following code in R (we highly recommend using R studio)

From the Command line
cd Folder_With_Tabulated_File
R
> library(ggplot2)
> Tab2<-read.table("Tabulated_File.txt", header=T)
> ggplot(Tab2, aes(Sample, Top1, ymin = Top2, ymax = Top1, fill = Mean_Top_Delta)) + geom_crossbar(stat='identity', linetype='blank', width=0.75) + xlab(label="Samples Tested") + ylab(label="Top Two Isoform Range") + theme(axis.text.x = element_text(angle=90,hjust = 0.9,vjust=0.5), axis.title = element_text(size=16, vjust=0.4)) + scale_y_continuous(limits=c(0,1), breaks=seq(0, 1, 0.1)) + scale_fill_gradient2("Mean\nInter\nIsoform\nDelta", low="blue", mid = "yellow", high="red", limits=c(0.0,1.0), midpoint = 0.5) + geom_hline(aes(yintercept=0.9), colour="blue", linetype="dashed")
> ggsave("PurityCheckerV2.png")
> q()

