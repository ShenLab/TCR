###### Boris Grinshpun 01/30/16 ###### 

# presets
rm(list=ls())
options(stringsAsFactors = FALSE)
setwd('~/Desktop/Sykes_alloresponse') # path to data files
outputpath="~/Desktop/Sykes_alloresponse/filtered/" # where filtered results are saved

### DATA LOADING ###

#read file
args <- commandArgs(TRUE)
file<-args[1]  # "HC10v12_count_nt.txt.unixed.txt"  # <- CHANGE THIS LINE TO FILENAME TO RUN WITHOUT COMMANDLINE INPUT
data=read.table(file,header=T)

# separate stim and unstim
N=names(data)
CD4names=grep('CD4',N)
if(length(CD4names)==3){
  CD4names=CD4names[c(1,3)]
}
CD8names=grep('CD8',N)
if(length(CD8names)==3){
  CD8names=CD8names[c(2,3)]
}
unstimnames=grep('unstim',N)
stimnames=grep('unstim',N,invert=T)

names(data)[intersect(CD4names,unstimnames)]<-"CD4.unstim"
names(data)[intersect(CD8names,unstimnames)]<-"CD8.unstim"
names(data)[intersect(CD4names,stimnames)]<-"CD4.stim"
names(data)[intersect(CD8names,stimnames)]<-"CD8.stim"
data=data[,c('nucleotide','total','CD4.unstim','CD8.unstim','CD4.stim','CD8.stim')]
data=data[rowSums(data[3:6],)>0,]
rownames(data)<-NULL
### END DATA LOADING ###


### FILTER FUNCTIONS ###
#flags
#0 - keep
#1 - discard (CD4 and CD8 same subtype unresolvable)
#2 - discard (CD4 and CD8 different subtype unresolvable -- weirdness)
#3 - refilter and keep (CD4 and CD8 same subtype unresolvable, but other subtype might shed light [combination of flags #1+#2] )
#4 - discard (no resolution at all)

sorting_filter<-function(X,cutoff=2){
  output=X
  
  Xfq=X[,c('CD4.unstim','CD8.unstim','CD4.stim','CD8.stim')]
  Xfq=sweep(Xfq,2,colSums(Xfq),'/') # normalize to frequencies
  ratios=t(apply(Xfq, 1, function(x) x/max(x))) # normalize so that highest frequency in each row is set to 1
  indices=data.frame(0+(ratios>=1/cutoff)) # index on whether or not the clone is below the frequency cutoff (0 if yes, 1 if no)
  indices$flag=rep(0,dim(indices)[1]) # initialize flags
 
  
  # collect indices for each subtype
  cd4stimind=which(indices$CD4.stim==1)
  cd4unstimind=which(indices$CD4.unstim==1)
  cd8stimind=which(indices$CD8.stim==1)
  cd8unstimind=which(indices$CD8.unstim==1)
  
  # determine flag which is then used to adjust indices for keeping or removing a clone.
  indices$flag[intersect(cd4stimind,cd8stimind)]=1 # same type (both stim)
  indices$flag[intersect(cd4unstimind,cd8unstimind)]=1  # same type (both unstim)
  indices$flag[intersect(cd4unstimind,cd8stimind)]=indices$flag[intersect(cd4unstimind,cd8stimind)]+2 # different type (one stim, one unstim, but different subset)
  indices$flag[intersect(cd8unstimind,cd4stimind)]=indices$flag[intersect(cd8unstimind,cd4stimind)]+2 # different type (one stim, one unstim, but different subset)
  properinds=which(indices$flag==0) # indices where to keep the clone
  
  # if there is only a single 1 in the row of indices because all other clones were below the cutoff threshhold, another 1 is added
  # to the same subset under the other condition. In other words, if only CD4stim is kept, add a 1 to allow CD4unstim.
  indices$CD4.unstim[intersect(properinds,cd4stimind)]=1
  indices$CD4.stim[intersect(properinds,cd4unstimind)]=1
  indices$CD8.unstim[intersect(properinds,cd8stimind)]=1
  indices$CD8.stim[intersect(properinds,cd8unstimind)]=1
  
  # discard indices based on flag as explained above.
  indices[indices$flag==1,1:4]=0
  indices[indices$flag==2,1:4]=0
  indices[indices$flag==4,1:4]=0
  
  refilter=which(indices$flag==3)
  for(i in refilter){
    if((indices[i,1]==1) & (indices[i,3]==1)){
      indices[i,2]=0
      indices[i,4]=0
      }
    if((indices[i,2]==1) & (indices[i,4]==1)){
      indices[i,1]=0
      indices[i,3]=0
      }
  }
  
  output[,3:6]=output[,3:6]*indices[,1:4]
  output$flag=indices$flag

  return(output)
  
}

reactive_filter<-function(X,cutoff=5){
      output=X
      
      Xfq=X[,c('CD4.unstim','CD8.unstim','CD4.stim','CD8.stim')]
      Xfq=sweep(Xfq,2,colSums(Xfq),'/') # normalize to frequencies
      
      Xfq$CD4keep=0+((Xfq$CD4.stim/Xfq$CD4.unstim)>=cutoff)
      Xfq$CD8keep=0+((Xfq$CD8.stim/Xfq$CD8.unstim)>=cutoff)
      Xfq[is.na(Xfq)]=0 # NA happens when the sequence is not present in both stim and unstim (e.g. present for CD4, but not in either CD8)
  
      output$CD4.stim=output$CD4.stim*Xfq$CD4keep
      output$CD8.stim=output$CD8.stim*Xfq$CD8keep
      return(output)
}

### END FILTER FUNCTIONS ###

#### APPLY FILTERING AND SAVE OUTPUT
data.f=sorting_filter(data)
data.ff=reactive_filter(data.f)

# recompute totals
data.ff$total=rowSums(data.ff[,c('CD4.unstim','CD8.unstim','CD4.stim','CD8.stim')])

# remove 0 rows
data.ff=data.ff[data.ff$total>0,]

write.table(data.ff,file=paste(outputpath,'/',file,'.filtered',sep=''),quote=F,sep='\t',row.names = F)
