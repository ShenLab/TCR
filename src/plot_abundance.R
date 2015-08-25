rm(list=ls())
setwd('/Users/bg2178/DropboxCGC/FarberSamples/productive')

donor<-'D201'
type<-'CD4'
tissue='naive'
#tissue<-"LLN-"

L<-Sys.glob(paste(donor,"*",sep=''))
L<-L[grep(type,L)]
L<-L[grep(tissue,L,ignore.case = T)]
L<-L[grep(69,L,invert=T)]

options(stringsAsFactors=F) 
require(plyr)

counts2hist<-function(c,id)
{
  c=c/sum(c)
  h<-table(c)
  h<-cbind(names(h),h)
  mode(h) <- "numeric"
  h<-as.data.frame.matrix(h)
  names(h)<-c("X","Y")
  h$id=as.character(id)
  return(h)
}


S=lapply(seq(L), function(i) counts2hist(read.table(L[i],header=T)$count,i))
S=ldply(S,rbind)

S1=subset(S,id==1)
S2=subset(S,id==2)
S3=subset(S,id==3)
#S4=subset(S,id==4)
#S5=subset(S,id==5)

png(file=paste('/Users/bg2178/DropboxCGC/newabundance/',donor,",",tissue,type,'.ycounts.png',sep=""),pointsize=5,res=1200,height=2000,width=2000)
plot(S1$X,S1$Y,xlim=c(min(S$X),max(S$X)),ylim=c(1,10^5),col='red',log='xy',main=paste(donor,", ",toupper(tissue),", ",type,sep=""),cex=1.1,
        cex.lab=1.2,cex.main=1.3,xlab='Copy number frequency',ylab='Number clones',asp=1,pch=20)

points(S2$X,S2$Y,col='blue',pch=20)
points(S3$X,S3$Y,col='green',pch=20)
#points(S4$X,S4$Y,col='lightblue',pch=20)
#points(S5$X,S5$Y,col='green',pch=20)


legend('topright',c("ILN","LLN","SP"),col=c('red','blue','green'),pch=c(20,20,20))
#legend('topright',c("69+","69-"),col=c('red','blue'),pch=c(20))
#legend('topright',c("ILN","LLN-1","LLN-2","LLN-3","SP"),col=c('red','blue','cyan','lightblue',"green"),pch=c(20))

dev.off()
