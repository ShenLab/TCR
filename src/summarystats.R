rm(list=ls())
options(stringsAsFactors=F) 
setwd('/Users/bg2178/DropboxCGC/FarberSamples/productive')

donor<-'D72'
type<-'CD4'
subset='TEM'

L<-Sys.glob(paste(donor,"*",sep=''))
L<-L[grep(type,L)]
L<-L[grep(subset,L,ignore.case = T)]
L<-L[grep(69,L,invert=T)]

S=lapply(seq(L), function(i) read.table(L[i],header=T))

calcH<-function(vals){
  fq=vals/sum(vals)
  H=-sum(fq*log2(fq))
  return(H)
}

calcC<-function(vals){
  H=calcH(vals)
  Hmax=log2(length(vals))
  C=1-H/Hmax
  return(C)
}

calcSI<-function(vals){
  fq=vals/sum(vals)
  si=sum(fq^2)
  return(si)
}

summarystats=data.frame(samples=L)
summarystats$nreads=sapply(S,function(x) sum(x$count))
summarystats$nclones=sapply(S,function(x) length(x$count))
summarystats$entropy=sapply(S, function(x) calcH(x$count))
summarystats$clonality=sapply(S, function(x) calcC(x$count))
summarystats$SI=sapply(S, function(x) calcSI(x$count))