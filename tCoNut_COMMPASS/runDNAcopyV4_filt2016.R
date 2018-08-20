#runDNAcopyExomeV2<-function(cghTSV){
  
library('DNAcopy')
Sys.setenv("DISPLAY"="hpc-hn01:0.0")

args = commandArgs(TRUE)
cghTSV = args[1]
pngName = args[2]  

cgh=read.table(cghTSV,header=TRUE,sep="\t")
  
fName=substr(cghTSV,1,nchar(cghTSV)-4)
CNA.object=CNA(cgh$Fold.Change,cgh$Chr,cgh$Position,data.type="logratio",sampleid=fName)

#segment.CNA.object=segment(CNA.object,alpha=0.001,verbose=1,undo.splits="sdundo",undo.SD=4,min.width=2)
segment.CNA.object=segment(CNA.object, alpha = 0.001, verbose = 1 , undo.splits = "sdundo", undo.SD = 3.8, min.width = 2 )
 
write.table(print(segment.CNA.object),file=paste(fName,'.seg',sep=''),sep="\t",row.names=FALSE)
  
png(file=paste(pngName,'.png',sep=''),width=500*4*3,height=500*2*3, res=300)
plot(segment.CNA.object,plot.type="w", ylim=c(-4,4))
dev.off()

for(chrom in 1:23)
#for(chrom in 1:22)
{
	png(file=paste(pngName,'.CHR.',chrom,'.png',sep=''),width=500*4*3,height=500*2*3, res=300)
	plot(subset(segment.CNA.object, chromlist=chrom),plot.type="c", ylim=c(-4,4))
	dev.off()
}

