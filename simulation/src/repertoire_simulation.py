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

def plot_abundance(x,y,a,k):
	'''
	x: copy number of largest clone
        y: number of unique clones with each copy number 1:x

	'''
	uniqueclones=y
	import matplotlib
	matplotlib.use('Agg')
	import matplotlib.pyplot as plt

	copyno=range(1,x+1)
	copyfq=[float(i)/sum(copyno) for i in copyno]
	
	fig=plt.figure()
	ax=plt.gca()
	eq="y=%.1e*x^%.2g" % (a,k)
	ax.scatter(copyfq,uniqueclones,lw = 0,label=eq)
	ax.set_yscale('log')
	ax.set_xscale('log')
	plt.axis([min(copyfq)*0.9,max(copyfq)*1.2,0.9, uniqueclones[0]*1.2])
	#plt.axis([0.9,(x+1)*1.2,0.9,uniqueclones[0]*1.2])
	plt.legend(loc='upper right',numpoints = 1)
	plt.savefig('testplot.png')
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
