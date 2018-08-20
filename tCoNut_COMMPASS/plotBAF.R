#plotBAF <- function(bafTXT,pngName){

Sys.setenv("DISPLAY"="hpc-hn01:0.0")

args <- commandArgs(TRUE)

bafTXT=args[1]
pngName=args[2]

#read in TSV files
baf<-read.table(bafTXT,header=TRUE,sep="\t")

baf[,2]=baf[,2]/1e6

#create PNG file
png(file=paste(pngName,'.png',sep=''),width=500*6*3,height=500*4*3, res=300)
par(mfrow=c(6,4))

for (i in unique(baf[,1])){
  tmpBAF <- subset(baf, Chromosome == i)
  
  plot(tmpBAF$Position,tmpBAF$BAF,
       type="p",
       col="dark blue",
       bg="dark blue",
       pch=16,
       ylim=c(0,1),
       yaxt="n",
       main=paste("Chromosome ",i),
       ylab="Frequencies",
       xlab="Physical Position (Mb)",
       xlim=c(floor(min(tmpBAF$Position)),ceiling(max(tmpBAF$Position))),
       xaxt="n",
       lwd=3.5,
       cex.lab=1.5,
       cex.main=1.5
  )
  abline(h=0.5,lty=2)
  
  #x=seq(floor(min(tmpCGH$Position)),ceiling(max(tmpCGH$Position)),25)
  y=seq(0,1,0.5)
  axis(2,at=y)
  x=seq(0,250,25)
  axis(1,at=x)
}
#xlim=c(floor(min(tmpCGH$Position)),ceiling(max(tmpCGH$Position))),
invisible(dev.off())

