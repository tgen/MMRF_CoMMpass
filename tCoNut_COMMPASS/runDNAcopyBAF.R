#runDNAcopyBAF<-function(bafTXT,fName){
 
Sys.setenv("DISPLAY"="hpc-hn01:0.0")

args <- commandArgs(TRUE)
bafTXT=args[1]
fName=args[2]
 
  library('DNAcopy')
  baf=read.table(bafTXT,header=TRUE,sep="\t")
  baf$BAF=abs(0.5-baf$BAF)
  
 ##fName=substr(bafTXT,1,nchar(bafTXT)-4)
  CNA.object=CNA(baf$BAF,baf$Chr,baf$Position,data.type="logratio",sampleid=fName)
  
  segment.CNA.object=segment(CNA.object,min.width=5)
  
  write.table(print(segment.CNA.object),file=paste(fName,'.seg',sep=''),sep="\t",row.names=FALSE)
  
  #pngName=substr(cghTSV,1,9)
  png(file=paste(fName,'.seg.png',sep=''),width=500*4*3,height=500*2*3, res=300)
  plot(segment.CNA.object,plot.type="w")
  dev.off()

