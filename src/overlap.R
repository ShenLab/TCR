rm(list=ls())
options(stringsAsFactors=F) 
setwd('/Users/bg2178/DropboxCGC/FarberSamples/productive')

# File specification
donor<-'D194'
type<-'CD4'
subset='naive'
top500=0
toppercent=1
percent=10^-4

L<-Sys.glob(paste(donor,"*",sep=''))
L<-L[grep(type,L)]
L<-L[grep(subset,L,ignore.case = T)]
L<-L[grep(69,L,invert=T)]

# Read the files into a data frame and combine replicates
require(plyr)

#Spleen
sp<-L[grep("SP",L)]
spfull=lapply(seq(sp),function(i) read.table(sp[i],header=T))
spfull=ldply(spfull,rbind)
Sf=data.frame(nucleotide=spfull$nucleotide,cts=spfull$count)
spfull=aggregate(cts ~ nucleotide, Sf, sum)
spfull=spfull[order(-spfull[,2]),]
spname=paste(unlist(strsplit(sp,'_'))[1:3],collapse='')
 # get top 500
if(top500==1 & length(spfull$cts)>=500 ){
ind=max(which(spfull$cts==(spfull$cts)[500]))
spfull=spfull[1:ind,]
}
if(toppercent==1){
  spfullfq=spfull$cts/sum(spfull$cts)
  ind=spfullfq>=percent
  spfull=spfull[ind,]
}


#LLN
lln<-L[grep("LLN",L)]
llnfull=lapply(seq(lln),function(i) read.table(lln[i],header=T))
llnfull=ldply(llnfull,rbind)
Lf=data.frame(nucleotide=llnfull$nucleotide,cts=llnfull$count)
llnfull=aggregate(cts ~ nucleotide, Lf, sum)
llnfull=llnfull[order(-llnfull[,2]),]
llnname=paste(unlist(strsplit(lln,'_'))[1:3],collapse='')      
# get top 500
if(top500==1 & length(llnfull$cts)>=500){
ind=max(which(llnfull$cts==(llnfull$cts)[500]))
llnfull=llnfull[1:ind,]
}
if(toppercent==1){
  llnfullfq=llnfull$cts/sum(llnfull$cts)
  ind=llnfullfq>=percent
  llnfull=llnfull[ind,]
}



#ILN
iln<-L[grep("ILN",L)]
ilnfull=lapply(seq(iln),function(i) read.table(iln[i],header=T))
ilnfull=ldply(ilnfull,rbind)
If=data.frame(nucleotide=ilnfull$nucleotide,cts=ilnfull$count)
ilnfull=aggregate(cts ~ nucleotide, If, sum)
ilnfull=ilnfull[order(-ilnfull[,2]),]
ilnname=paste(unlist(strsplit(iln,'_'))[1:3],collapse='')      
# get top 500
if(top500==1 & length(ilnfull$cts)>=500 ){
ind=max(which(ilnfull$cts==(ilnfull$cts)[500]))
ilnfull=ilnfull[1:ind,]
}
if(toppercent==1){
  ilnfq=ilnfull$cts/sum(ilnfull$cts)
  ind=ilnfq>=percent
  ilnfull=ilnfull[ind,]
}


R<-list()
R[[1]]=ilnfull
R[[2]]=llnfull
R[[3]]=spfull
# Merge the data frames
Q=R[[1]]
for(i in 1:(length(R)-1)){
  Q<-merge(Q,R[[i+1]],by='nucleotide',all=T)
}
Q[is.na(Q)]=0
names(Q)[2]<-ilnname
t1=unlist(strsplit(iln,'_'))[2]
names(Q)[3]<-llnname
t2=unlist(strsplit(lln,'_'))[2]
names(Q)[4]<-spname
t3=unlist(strsplit(sp,'_'))[2]

## Get the overlap
uq1=sum(0+(Q[2]>0 & Q[3]==0 & Q[4]==0))
uq2=sum(0+(Q[3]>0 & Q[2]==0 & Q[4]==0))
uq3=sum(0+(Q[4]>0 & Q[3]==0 & Q[2]==0))
o1=sum(0+(Q[2]>0 & Q[3]>0 & Q[4]==0))
o2=sum(0+(Q[2]>0 & Q[4]>0 & Q[3]==0))
o3=sum(0+(Q[3]>0 & Q[4]>0 & Q[2]==0))
otot=sum(0+(Q[2]>0 & Q[3]>0 & Q[4]>0))

print(paste('Unique ',names(Q)[2],' = ',uq1,sep=''))
print(paste('Unique ',names(Q)[3],' = ',uq2,sep=''))
print(paste('Unique ',names(Q)[4],' = ',uq3,sep=''))
print(paste('Overlap ',names(Q)[2],'-', names(Q)[3],' = ',o1,sep=''))
print(paste('Overlap ',names(Q)[2],'-', names(Q)[4],' = ',o2,sep=''))
print(paste('Overlap ',names(Q)[3],'-', names(Q)[4],' = ',o3,sep=''))
print(paste('Overlap ALL = ',otot,sep=''))


# Inverse Jaccard Ind
overlap=(uq1+uq2+uq3)/(uq1+uq2+uq3+o1+o2+o3+otot)
print(overlap)
print(uq1+uq2+uq3+o1+o2+o3+otot)


# make venn diagrams
require(VennDiagram)

if (top500==1){
  venn.diagram(list(R[[1]]$nucleotide,R[[2]]$nucleotide,R[[3]]$nucleotide),fill = c("red", "blue","green"),
               alpha = c(0.5, 0.5,0.5), category.names = c(t1,t2,t3),cat.fontface = 4,lty=2, imagetype = "png",main=paste(donor,", ",subset,", ",type,sep=""),main.fontface=4,resolution = 750,height = 3000, width = 3000,asp=1,cex=1.1,main.cex=1.4,filename=paste("/Users/bg2178/DropboxCGC/overlap_venn/",donor,",",subset,",",type,"_top500.png",sep=""))
} else {
  venn.diagram(list(R[[1]]$nucleotide,R[[2]]$nucleotide,R[[3]]$nucleotide),fill = c("red", "blue","green"),
             alpha = c(0.5, 0.5,0.5), category.names = c(t1,t2,t3),cat.fontface = 4,lty=2, imagetype = "png",main=paste(donor,", ",subset,", ",type,sep=""),main.fontface=4,resolution = 750,height = 3000, width = 3000,asp=1,cex=1.1,main.cex=1.4,filename=paste("/Users/bg2178/DropboxCGC/overlap_venn/",donor,",",subset,",",type,".png",sep=""))
}
