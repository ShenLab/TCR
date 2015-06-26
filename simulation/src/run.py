from repertoire_simulation import make_pwlaw, make_repertoire
from diversitymeasures import calc_entropy as entropy
from diversitymeasures import make_hist
import itertools

k=-2.5
a=1e7
s=1e8

x,y=make_pwlaw(k,a,s)
repertoire=make_repertoire(x,y,650)
del y

# x,y expanded
expandedclones=[[i]*len(repertoire[i]) for i in repertoire.keys()]
expandedclones=list(itertools.chain(*expandedclones))
allVJ=list(itertools.chain(*repertoire.values()))

vjcounts=zip(*make_hist(allVJ))

Hcdr3=entropy(expandedclones)
Hvj=entropy(vjcounts[1])

print Hcdr3
print Hvj

