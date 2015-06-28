#!/usr/bin/env python

def make_pwlaw(slope,intercept,nseqs):
	'''
	power law: y=intercept*x^slope
	nseqs : approximate total number of cells

	'''
	from math import floor
	n=0
	x=1
	y=[]
	while n<nseqs:
		yy=floor(intercept*x**slope)
		if yy==0:
			break;
		n=n+yy*x
		y.append(int(yy))
		x=x+1
	return x-1,y
		
def make_repertoire(x,y,usage):
	'''
	x: copy number of largest clone
	y: number of unique clones with each copy number 1:x
	usage: number of possible VJ combinations
	
	'''
	import numpy as np
	from random import randint

	copyno=range(1,x+1)
	VJlabel=[[randint(1, usage) for i in range(j)] for j in y]
	repertoire=dict(zip(copyno,VJlabel))
	return repertoire

def sample_repertoire(repertoire, samplesize):
	'''
	sample from a repertoire
	'''
	from random import sample
	import itertools
	from diversitymeasures import make_hist

	# expand the repertoire and index clones
	expandedVJ=[[(n,j) for j in vj] for n,vj in repertoire.items()]
	expandedVJ=list(itertools.chain(*expandedVJ))
	expandedVJ=[(i,n,vj) for i,(n,vj) in enumerate(expandedVJ)]  #i==cloneid,n==copyno,j==VJid
	
	fullCDR3=[[[i,vj]]*n for i,n,vj in expandedVJ]
	fullCDR3=tuple(itertools.chain(*fullCDR3))
	sampleCDR3=sample(fullCDR3,samplesize)
	
	h=make_hist(sampleCDR3)

	newrepertoire=dict()
	# convert into repertoire format [copy # : VJ list]
	for (i,j),k in h:
		newrepertoire.setdefault(k,[]).append(j)

	return newrepertoire


def sample_nb(repertoire,p,cov):
	from numpy.random import negative_binomial as nb
			
	sample=[[(nb(i*cov,p),j[k]) for k in range(len(j))] for i,j in repertoire.items()]
	return sample
		

def plot_abundance(x,y,leg,name='output',ymax=None,col='b'):
	'''
	x: copy number of largest clone
        y: number of unique clones with each copy number 1:x

	'''
	uniqueclones=y
	import matplotlib
	matplotlib.use('Agg')
	import matplotlib.pyplot as plt

	copyno=x
	copyfq=[float(i)/sum(copyno) for i in copyno]
	fig=plt.figure()
	ax=plt.gca()
	#eq="y=%.1e*x^%.2g" % (a,k)
	ax.scatter(copyfq,uniqueclones,c=col,lw = 0,label=leg)
	ax.set_yscale('log')
	ax.set_xscale('log')

	if ymax is None:
		ymax=max(uniqueclones)*1.2	
	plt.axis([min(copyfq)*0.9,max(copyfq)*1.2,0.9, ymax])
	#plt.axis([0.9,(x+1)*1.2,0.9,uniqueclones[0]*1.2])
	plt.legend(loc='upper right',numpoints = 1)
	plt.savefig('%s.png' % name )
	# plt.show()

def plot_abundance_pair(x1,y1,x2,y2,leg=None,name='output',ymax=None):
        '''
        x: copy number of largest clone
        y: number of unique clones with each copy number 1:x

        '''
	cellcopy=x1
        uniquecells=y1
	readcopy=x2
	uniqueclones=y2
	
	cellfq=[float(i)/sum(cellcopy) for i in cellcopy]
	readfq=[float(i)/sum(readcopy) for i in readcopy]

        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
	from copy import deepcopy

        fig=plt.figure()
	ax=fig.add_subplot(111)
	ax.scatter(cellfq, uniquecells, c='r',marker='^',label='cells',lw=0)
	ax.scatter(readfq, uniqueclones, c='b', marker='o',label='reads',lw = 0)
        #ax.scatter(copyfq,uniqueclones,lw = 0,label=leg)
        ax.set_yscale('log')
        ax.set_xscale('log')
        if ymax is None:
                ymax=max(uniquecells)*1.2
	
	allfq=deepcopy(cellfq)
	allfq.extend(readfq)
        plt.axis([min(allfq)*0.9,max(allfq)*1.2,0.9, ymax])
        #plt.axis([0.9,(x+1)*1.2,0.9,uniqueclones[0]*1.2])
        plt.legend(loc='upper right',numpoints = 1)
        plt.savefig('%s.png' % name )
        # plt.show()


def main():
        import os
        from argparse import ArgumentParser # for specifying input arguments

        def get_args():
                '''
                Collect input arguments

                '''
                # Usage:
                parser=ArgumentParser(description="Generate power law distribution", usage="%(prog)s  sim_pwlaw.py [-h] [-k slope] [-a intercept] [-v number of VJ pairs]" )

                # Positional arguments
                parser.add_argument("-k", help="power law slope (y=ax^k)")
                parser.add_argument("-a" ,help="power law intercept (y=ax^k)")
                parser.add_argument("-v", default=650, help="Number of allowed VJ combinations")
		parser.add_argument("-s", default=1e9, help="number of sequences")

                return parser

        parser=get_args() # get input arguments
        args=parser.parse_args()
 
        if args.a is None:
                print "MISSING INTERCEPT -a"
                exit()
	else: a=float(args.a)
        if args.k is None:
                print "MISSING SLOPE -k"
                exit()
	else: k=float(args.k) 
	vj=int(args.v)
	s=float(args.s)

	#x,y=make_pwlaw(k,a,s) # construct power law

	#repertoire=make_repertoire(x,y,vj) # build clonal repertoire	
	#plot_abundance(x,y,a,k)
		

if __name__ == '__main__':
	''' 
	python sim_pwlaw.py -h to learn more

	'''
	main()
