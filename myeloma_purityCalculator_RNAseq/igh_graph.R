#######################################################################
#######################################################################
## RNAseq based Myeloma Purity Calculator
##
## Copyright (c) 2017 Translational Genomics Research Institute
##
## This software may be modified and distributed under the terms
## of the MIT license.  See the LICENSE file for details.
##
## Major Contributors: Jonathan J. Keats 
## Minor Contributors:
#######################################################################
#######################################################################

#Load the ggplot library
#library(ggplot2)
if (!require('ggplot2')) install.packages('ggplot2'); library('ggplot2')

#Read in the IgH table
IGH<-read.table("Graph_IgH.txt", header=T)
TITLE<-readChar("title.txt",nchars=1e6)
head(TITLE)

#levels=c("IGHM","IGHD","IGHG3","IGHG1","IGHA1","IGHG2","IGHG4","IGHE","IGHA2","IGKC","IGLC1","IGLC2","IGLC3","IGLC7")

#Generate IgH graph
ggplot(data=IGH, aes(x = CommonName, y=Percentage, fill = TotalFrequency)) + 
	geom_bar(stat='identity', width=.5) + 
	xlab(label="Immunoglobulin Heavy Chain Isoforms") + 
	ylab(label="Respective IgH Read Fraction") + 
	theme(axis.text = element_text(size=10), axis.text.x = element_text(size=9, angle=90,hjust = 0.9, vjust = 0.5), axis.title = element_text(size=14, vjust=0.4)) + 
	scale_y_continuous(limits=c(0,1), breaks=seq(0, 1, 0.1)) + 
	scale_fill_gradient2("Percent\nTotal\nReads", low="blue", mid = "yellow", high="red", limits=c(0.0,0.1), midpoint = 0.05, na.value = "red") + 
	geom_hline(aes(yintercept=0.9), colour="Red", linetype="dashed") + 
#	annotate("pointrange", x = "IGHM", y = 0.4, ymin = 0.3, ymax = 0.5, color = "red", size = 0.5) + 
	facet_grid(. ~ Locus, scales="free", space="free") +
	ggtitle(TITLE)
ggsave(file="IGH.png", width = 10, height = 6)


#Read in the IgL Table
IGL<-read.table("Graph_IgL.txt", header=T)

#Generate IgL Variable graph
ggplot(data=IGL, aes(x = CommonName, y = Percentage, fill = TotalFrequency)) + 
	geom_bar(stat='identity', width=.5) + 
	xlab(label="Immunoglobulin Light Chain Isoforms") + 
	ylab(label="Respective IgL Read Fraction") + 
	theme(axis.text = element_text(size=10), axis.text.x = element_text(size=7, angle=90,hjust = 0.9, vjust = 0.5), axis.title = element_text(size=14, vjust=0.4)) + 
	scale_y_continuous(limits=c(0,1), breaks=seq(0, 1, 0.1)) + 
	scale_fill_gradient2("Percent\nTotal\nReads", low="blue", mid = "yellow", high="red", limits=c(0.0,0.1), midpoint = 0.05, na.value = "red") + 
	geom_hline(aes(yintercept=0.9), colour="Red", linetype="dashed") + 
#	coord_flip() + 
	facet_grid(~ Locus, scales="free", space="free") + 
	ggtitle(TITLE)

ggsave(file="IGL.png", width = 10, height = 6)

#Quit R Session
q()
