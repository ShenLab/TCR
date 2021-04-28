from repertoire_simulation import * 
from diversitymeasures import calc_entropy as entropy
from diversitymeasures import make_hist
import itertools

k=-1.5
a=20000
s=700000
mm=10**5
x,y=make_pwlaw(k,a,s)
eq="y=%.1e*x^%.2g" % (a,k)
repertoire=make_repertoire(x,y,650)
x=xrange(1,x+1)
#plot_abundance(range(1,x+1),y,eq,'fullrep',ymax=mm,col='r')
print sum([i*j for i,j in zip(x,y)])

frac_clones=0.0001
cells_per_clone=0.7
s=0.2
expandedrepertoire=expand_repertoire(repertoire,frac_clones,cells_per_clone,s)

if 0 in expandedrepertoire.keys():
        expandedrepertoire.pop(0)

y=[len(j) for j in expandedrepertoire.values()]
x=expandedrepertoire.keys()

plot_abundance(x,y,eq,'expanded',ymax=mm,col='r')
print sum([i*j for i,j in zip(x,y)])

sample=sample_nb(expandedrepertoire,0.6,5)
sample=list(itertools.chain(*sample))

# convert sample into repertoire format <-- added to simulation code
newrepertoire=dict()
for i,j in sample:
	newrepertoire.setdefault(i,[]).append(j)

if 0 in newrepertoire.keys():
	newrepertoire.pop(0)

y2=[len(j) for j in newrepertoire.values()]
x2=newrepertoire.keys()

#plot_abundance(x2,y2,'nb','negative_bin',ymax=mm)
#x=range(1,x+1)
plot_abundance_pair(x,y,x2,y2,name='plotpair_final',ymax=mm)

print sum([i*j for i,j in zip(x2,y2)])
print sum(y2)
print "\n"


exit()
# Convert to diversity format
def expand(rep):
	expandedclones=[[i]*len(rep[i]) for i in rep.keys()]
	expandedclones=list(itertools.chain(*expandedclones))

	# i == clone id, n==copy no, j==VJid
	expandedrep=[[(n,j) for j in vj] for n,vj in rep.items()]
	expandedrep=list(itertools.chain(*expandedrep))
	expandedrep=[(i,j,k) for i,(j,k) in enumerate(expandedrep)]
	fullCDR3=[[[i,vj]]*n for i,n,vj in expandedrep]
        expandedVJ=zip(*tuple(itertools.chain(*fullCDR3)))[1]
	expandedVJ=zip(*make_hist(expandedVJ))[1]

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
