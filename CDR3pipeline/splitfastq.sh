#!/bin/bash

# Take a sample and split it up into piece of size 2000000 for quicker processing

sample=$1 # sample name
indir=$2 # sample directory
outdir=$3 # directory to place the output

mkdir -p $outdir  # make output directory if it doesn't exist

# split into pieces
echo "split -l 10000000 $indir/$sample $outdir/`basename $sample .fastq`.part." | qsub -N $sample -cwd -l mem=2G,time=3:: -o log -e log
