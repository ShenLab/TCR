#!/bin/bash

prefix=$1 # file prefix (followed by aa,ab,ac parts)
datadir=$2 # directory where split files are located
fastaref=$3 # path to the reference fasta file
outdir=$4 # output directory for bams

# current human fasta can be found here: /ifs/scratch/c2b2/ys_lab/yshen/TCR/GBM_LGG/fastas/human_g1k_v37.MaskTRB.TRBp8.fasta
 

mkdir -p $outdir

fqfiles=$datadir/$prefix*.part*
for fq in $fqfiles
do
	qsub -N bwa`basename $fq` -cwd -pe smp 4 -R y -l mem=5G,time=2:: bwa_job.sh $fq $fastaref $outdir

done
