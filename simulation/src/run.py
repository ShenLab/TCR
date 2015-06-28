from repertoire_simulation import * 
from diversitymeasures import calc_entropy as entropy
from diversitymeasures import make_hist
import itertools

k=-2
a=50
s=20000

x,y=make_pwlaw(k,a,s)
eq="y=%.1e*x^%.2g" % (a,k)
repertoire=make_repertoire(x,y,650)
plot_abundance(range(1,x+1),y,eq,'fullrep2')

del y

sample=sample_nb(repertoire,0.8,20)
sample=list(itertools.chain(*sample))

newrepertoire=dict()
for i,j in sample:
	newrepertoire.setdefault(i,[]).append(j)

newrepertoire.pop(0)

y=[len(j) for j in newrepertoire.values()]
x=newrepertoire.keys()

print sum([i*j for i,j in zip(x,y)])

#sample=list(itertools.chain(*sample))
#print sample

#s=make_hist(sample)
#print s

plot_abundance(x,y,'nb','negative_bin')

#print sum([i*j for i,j in zip(*s)])

#sample=sample_repertoire(repertoire,int(1e6))

#x=sample.keys()
#y=[len(j) for i,j in sample.items()]
#plot_abundance(s[0],s[1],'nb','negative_bin0.6')

# x,y expanded
#expandedclones=[[i]*len(repertoire[i]) for i in repertoire.keys()]
#expandedclones=list(itertools.chain(*expandedclones))

# i == clone id, n==copy no, j==VJid

#expandedVJ=[[(n,j) for j in vj] for n,vj in repertoire.items()]
#expandedVJ=list(itertools.chain(*expandedVJ))
#expandedVJ=[(i,j,k) for i,(j,k) in enumerate(expandedVJ)]
#print expandedVJ[0:5]
#print len(expandedVJ)

#allVJ=list(itertools.chain(*repertoire.values()))
#expandedVJ=[[j]*i for i,j in expandedVJ]
#allVJ=list(itertools.chain(*expandedVJ))
#print expandedVJ[0:5]


#Hcdr3=entropy(expandedclones)
#Hvj=entropy(vjcounts[1])

#print Hcdr3
#print Hvj

