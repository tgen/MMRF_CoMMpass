This is designed to determine if an RNAseq result bam file contains a pure monoclonal B cell population

REQUIREMENTS

1) samtools
2) featurecounts (subread package)
3) R
4) ggplot2 library in R

All other required files for GRCh37 are provided but you can build the required annotation files by changing the BUILD_FILES flag to "YES"

Currently, the application will process all .bam files in a folder from which the script is called

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


TEST
