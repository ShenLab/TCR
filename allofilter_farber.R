###### Boris Grinshpun 01/30/16 ###### 

# presets
rm(list=ls())
options(stringsAsFactors = FALSE)
setwd('~/DropboxCGC/FarberSamples/productive') # path to data files
outputpath="~/DropboxCGC/FarberSamples/productive/filtered/" # where filtered results are saved

### DATA LOADING ###

#read file
args <- commandArgs(TRUE)
donor=args[1]
files=Sys.glob(paste(donor,'*.tsv',sep=''))
data=lapply(files, function(f) read.table(f,header=T))
data=lapply(data,function(x) {rownames(x)=x$nucleotide;x;})
mynames=sapply(files, function(f) strsplit(f,split='.',fixed=T)[[1]][1])

data.simplified=lapply(seq(data), function(i){x=data.frame(data[[i]]$nucleotide,data[[i]]$count); colnames(x)=c("nucleotide",mynames[i]);x;})
dd=Reduce(function(...) merge(..., all=T,by='nucleotide'), data.simplified)
dd[is.na(dd)]=0


sorting_filter<-function(X,cutoff=2){
  output=data.frame(nucleotide=X$nucleotide)
  
  Xfq=X[,2:length(X)]
  Xfq=sweep(Xfq,2,colSums(Xfq),'/') # normalize to frequencies
  ratios=t(apply(Xfq, 1, function(x) x/max(x))) # normalize so that highest frequency in each row is set to 1
  indices=data.frame(0+(ratios>=1/cutoff)) # index on whether or not the clone is below the frequency cutoff (0 if yes, 1 if no)
 # indices$flag=rep(0,dim(indices)[1]) # initialize flags
  
  CD4vals=data.frame(CD4=indices[,grep('CD4',names(indices),ignore.case=T,value=T)])
  CD8vals=data.frame(CD8=indices[,grep('CD8',names(indices),ignore.case=T,value=T)])
  CD4ratios=data.frame(CD4=ratios[,grep('CD4',colnames(ratios),ignore.case=T,value=T)])
  CD8ratios=data.frame(CD8=ratios[,grep('CD8',colnames(ratios),ignore.case=T,value=T)])  
  
  if(dim(CD4vals)[2]>1 & dim(CD8vals)[2]>1){
  CD4vals$total=rowSums(CD4vals)
  CD4ratios$total=rowSums(0+(CD4ratios>0))
  CD8vals$total=rowSums(CD8vals)
  CD8ratios$total=rowSums(0+(CD8ratios>0))
  } else {CD4vals$total=CD4vals$CD4;CD4ratios$total=0+(CD4ratios>0); CD8vals$total=CD8vals$CD8; CD8ratios$total=0+(CD8ratios>0)}
  
  indices$CD4present=0+(CD4vals$total>0)
  indices$CD8present=0+(CD8vals$total>0)
  
  # resolve ambiguity to the best possible first by indices, then by ratios across all samples.
  ambiguous=which(indices$CD4present>0 & indices$CD8present>0)
  indices$CD4present[ambiguous[CD4vals$total[ambiguous]/CD8vals$total[ambiguous]<=0.5]]=0
  indices$CD8present[ambiguous[CD4vals$total[ambiguous]/CD8vals$total[ambiguous]>=2]]=0
  
  ambiguous=which(indices$CD4present>0 & indices$CD8present>0)
  indices$CD4present[ambiguous[CD4ratios$total[ambiguous]/CD8ratios$total[ambiguous]<=0.5]]=0
  indices$CD8present[ambiguous[CD4ratios$total[ambiguous]/CD8ratios$total[ambiguous]>=2]]=0
  
  ambiguous=which(indices$CD4present>0 & indices$CD8present>0)
  indices$CD4present[ambiguous]=0
  indices$CD8present[ambiguous]=0
  
  output$CD4inds=indices$CD4present
  output$CD8inds=indices$CD8present
  
  return(output)
  
}

inds=sorting_filter(dd)
rownames(inds)<-inds$nucleotide
data.f=data

# filter
for(i in seq(data)){
    n=mynames[i]
    d=data[[i]]
    seq=intersect(rownames(d),rownames(inds))
    d=d[seq,]
    currentinds=inds[seq,]
    if(length(grep("CD4",n))){
      d$count= d$count*currentinds$CD4inds
    } else if(length(grep("CD8",n))){
      d$count= d$count*currentinds$CD8inds
    } else{print('PROBLEM')}
    data.f[[i]]=d
}

# remove 0s
data.ff=lapply(data.f, function(x) x[x$count>0,])

for(i in seq(data.ff)){
  write.table(data.ff[[i]],file=paste(outputpath,'/',mynames[i],'.filtered',sep=''),quote=F,sep='\t',row.names = F)
}
