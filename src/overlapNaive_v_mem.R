rm(list=ls())
options(stringsAsFactors=F) 
setwd('/Users/bg2178/DropboxCGC/FarberSamples/productive')

#require(sm)

donor<-'D201'
type<-'CD8'
tissue1='SP'
tissue2='ILN'
incr=10

L<-Sys.glob(paste(donor,"*",sep=''))
L<-L[grep(type,L)]
Lmem<-L[grep("tem",L,ignore.case=T)]
Lnaive<-L[grep("naive",L,ignore.case=T)]

# GET OVERLAP FOR MEM
fmem1=read.table(Lmem[grep(tissue1,Lmem,ignore.case=T)],header=T)
fmem2=read.table(Lmem[grep(tissue2,Lmem,ignore.case=T)],header=T)
fmem1$fq=fmem1$count/sum(fmem1$count)
fmem2$fq=fmem2$count/sum(fmem2$count)

# separate by count
mem1uniqcount=unique(fmem1$count)
mem1clones=lapply(mem1uniqcount, function(x)  fmem1$nucleotide[fmem1$count==x] )

mem2uniqcount=unique(fmem2$count)
mem2clones=lapply(mem2uniqcount, function(x)  fmem2$nucleotide[fmem2$count==x] )

# get overlap for mem
r1mem=unlist(lapply(seq(mem1clones), function(i) length(intersect(mem1clones[[i]],fmem2$nucleotide))))
r2mem=unlist(lapply(seq(mem2clones), function(i) length(intersect(mem2clones[[i]],fmem1$nucleotide))))

r11mem=sapply(seq(from=1, to=length(r1mem),by=(incr)), function(i) sum(r1mem[i:(i+incr)])/length(unlist(mem1clones[i:(i+incr)])))
mem11counts=sapply(seq(from=1, to=length(r1mem),by=(incr)),function(i) mean(mem1uniqcount[i:(i+incr)]))
mem11fq=mem11counts/sum(fmem1$count)
mem11fq=mem11fq[!is.na(mem11fq)]
r11mem=r11mem[!is.na(r11mem)]

sapply(seq(from=1, to=length(r2mem),by=(incr)), function(i) c(i, i+incr))
r22mem=sapply(seq(from=1, to=length(r2mem),by=(incr)), function(i) sum(r2mem[i:(i+incr)])/length(unlist(mem2clones[i:(i+incr)])))
mem22counts=sapply(seq(from=1, to=length(r2mem),by=(incr)),function(i) mean(mem2uniqcount[i:(i+incr)]))
mem22fq=mem22counts/sum(fmem2$count)
mem22fq=mem22fq[!is.na(mem22fq)]
r22mem=r22mem[!is.na(r22mem)]

# separate by frequency
#mem1uniqfq=unique(fmem1$fq)
#mem1clones=lapply(mem1uniqfq, function(x)  fmem1$nucleotide[fmem1$fq==x] )

#mem2uniqfq=unique(fmem2$fq)
#mem2clones=lapply(mem2uniqfq, function(x)  fmem2$nucleotide[fmem2$fq==x] )

# get overlap for mem
#r1mem=unlist(lapply(seq(mem1clones), function(i) (length(intersect(mem1clones[[i]],fmem2$nucleotide))/length(mem1clones[[i]]))))
#r2mem=unlist(lapply(seq(mem2clones), function(i) length(intersect(mem2clones[[i]],fmem1$nucleotide))/length(mem2clones[[i]])))



#plot(seqrange1,r1mem,log='x',pch=20)


# GET OVERLAP FOR NAIVE
fnaive1=read.table(Lnaive[grep(tissue1,Lnaive,ignore.case=T)][1],header=T)
fnaive2=read.table(Lnaive[grep(tissue2,Lnaive,ignore.case=T)][1],header=T)

# separate by count
naive1uniqcount=unique(fnaive1$count)
naive1clones=lapply(naive1uniqcount, function(x)  fnaive1$nucleotide[fnaive1$count==x] )

naive2uniqcount=unique(fnaive2$count)
naive2clones=lapply(naive2uniqcount, function(x)  fnaive2$nucleotide[fnaive2$count==x] )

# get overlap for naive
r1naive=unlist(lapply(seq(naive1clones), function(i) length(intersect(naive1clones[[i]],fnaive2$nucleotide))))
r2naive=unlist(lapply(seq(naive2clones), function(i) length(intersect(naive2clones[[i]],fnaive1$nucleotide))))

sapply(seq(from=1, to=length(r1naive),by=(incr)), function(i) c(i, i+incr))
r11naive=sapply(seq(from=1, to=length(r1naive),by=(incr+1)), function(i) sum(r1naive[i:(i+incr)])/length(unlist(naive1clones[i:(i+incr)])))
naive11counts=sapply(seq(from=1, to=length(r1naive),by=(incr+1)),function(i) mean(naive1uniqcount[i:(i+incr)]))
naive11fq=naive11counts/sum(fnaive1$count)
naive11fq=naive11fq[!is.na(naive11fq)]
r11naive=r11naive[!is.na(r11naive)]

r22naive=sapply(seq(from=1, to=length(r2naive),by=(incr+1)), function(i) sum(r2naive[i:(i+incr)])/length(unlist(naive2clones[i:(i+incr)])))
naive22counts=sapply(seq(from=1, to=length(r2naive),by=(incr+1)),function(i) mean(naive2uniqcount[i:(i+incr)]))
naive22fq=naive22counts/sum(fnaive2$count)
naive22fq=naive22fq[!is.na(naive22fq)]
r22naive=r22naive[!is.na(r22naive)]

#fnaive1$fq=fnaive1$count/sum(fnaive1$count)
#fnaive2$fq=fnaive2$count/sum(fnaive2$count)

# separate by frequency
#naive1uniqfq=unique(fnaive1$fq)
#naive1clones=lapply(naive1uniqfq, function(x)  fnaive1$nucleotide[fnaive1$fq==x] )
#naive2uniqfq=unique(fnaive2$fq)
#naive2clones=lapply(naive2uniqfq, function(x)  fnaive2$nucleotide[fnaive2$fq==x] )

# get overlap for naive
#r1naive=unlist(lapply(seq(naive1uniqfq), function(i) length(intersect(naive1clones[[i]],fnaive2$nucleotide))/length(naive1clones[[i]])))
#r2naive=unlist(lapply(seq(naive2uniqfq), function(i) length(intersect(naive2clones[[i]],fnaive1$nucleotide))/length(naive2clones[[i]])))



alloverlap=c(r11mem,r22mem,r11naive,r22naive)
#allfq=c(mem1uniqfq,mem2uniqfq,naive1uniqfq,naive2uniqfq)
allfq=c(mem11fq,mem22fq,naive11fq,naive22fq)

# PLOT
#png(paste('/Users/bgrinshpun/DropboxCGC/circosgo/scatteroutput/',donor,'_',tissue1,'_',type,tissue2,'.png',sep=''))
#par(bg=NA) 

plot(mem11fq,r11mem,ylim=c(0,max(alloverlap)),xlim=c(min(allfq),max(allfq)*1.2),pch=20,cex.axis=1.4,cex=1.6,cex.lab=1.5,xlab="clone frequency",ylab="overlap frequency",log='x')
points(mem22fq,r22mem,pch=20,col='blue',cex=1.5)
points(naive11fq,r11naive,pch=20,col='red',cex=1.5)
points(naive22fq,r22naive,pch=20,col='green',cex=1.5)


#plot(mem1uniqfq,r1mem,ylim=c(0,max(alloverlap)),xlim=c(min(allfq),max(allfq)*1.2),pch=20,xlab="clone frequency",ylab="overlap frequency",log='x')
#points(mem2uniqfq,r2mem,col='blue',pch=20)
#points(naive1uniqfq,r1naive,col='red',pch=18)
#points(naive2uniqfq,r2naive,col='green',pch=18)
legend('topleft',c(paste(tissue1,'mem',paste=''),
                    paste(tissue2,'mem',paste=''),paste(tissue1,'naive',paste=''),paste(tissue2,'naive',paste='')),pch=20,col=c('black','blue','red','green'),cex=1.3, pt.cex=1.3)
#dev.off()


