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


### Concept

We expect that a pure myeloma cell population should have a monoclonal expression of a single heavy and light chain molecule. We assume the unique heavy chain variable and constant gene segments will represent the vast majority of all RNAseq reads aligned to the IgH locus and similarly as single light chain variable element and constant region will represent the majority of the light chain alignments.  We expect that a heavy chain and light chain expressing patient will have single IgH VH and CH high frequency alleles and single high frequencyVK or VL and CK or CL alleles as shown.


###### Expected IgH Pattern
![IgH1](https://github.com/tgen/MMRF_CoMMpass/tree/master/myeloma_purityCalculator_RNAseq/images/Expected_Good_Myeloma_IgH.png)

testing
![IgH Test](/myeloma_purityCalculator_RNAseq/images/Typical_IgH.jpg)

###### Expected IgL Pattern
![IgL1](/images/Expected_Good_Myeloma_IgL.png)

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

