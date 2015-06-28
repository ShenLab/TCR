from repertoire_simulation import * 
from diversitymeasures import calc_entropy as entropy
from diversitymeasures import make_hist
import itertools

k=-1.5
a=500
s=200000
mm=10**4
x,y=make_pwlaw(k,a,s)
eq="y=%.1e*x^%.2g" % (a,k)
repertoire=make_repertoire(x,y,650)
plot_abundance(range(1,x+1),y,eq,'fullrep',ymax=mm,col='r')
print sum([i*j for i,j in zip(range(1,x+1),y)])
#del y

sample=sample_nb(repertoire,0.6,10)
sample=list(itertools.chain(*sample))

# convert sample into repertoire format <-- added to simulation code
newrepertoire=dict()
for i,j in sample:
	newrepertoire.setdefault(i,[]).append(j)

if 0 in newrepertoire.keys():
	newrepertoire.pop(0)

y2=[len(j) for j in newrepertoire.values()]
x2=newrepertoire.keys()

plot_abundance(x2,y2,'nb','negative_bin',ymax=mm)
#x=range(1,x+1)
plot_abundance_pair(x,y,x2,y2,name='plotpair',ymax=mm)

print sum([i*j for i,j in zip(x2,y2)])
print sum(y2)

# Convert to diversity format
def expand(rep):
	expandedclones=[[i]*len(rep[i]) for i in rep.keys()]
	expandedclones=list(itertools.chain(*expandedclones))

	# i == clone id, n==copy no, j==VJid
	expandedVJ=[[(n,j) for j in vj] for n,vj in rep.items()]
	expandedVJ=list(itertools.chain(*expandedVJ))
	expandedVJ=[(i,j,k) for i,(j,k) in enumerate(expandedVJ)]

	return expandedclones, expandedVJ


exC,exVJ=expand(repertoire)
Hcdr3=entropy(exC)
Hvj=entropy(exVJ)

print Hcdr3
print Hvj


exC,exVJ=expand(newrepertoire)
Hcdr3=entropy(exC)
Hvj=entropy(exVJ)

print Hcdr3
print Hvj 
