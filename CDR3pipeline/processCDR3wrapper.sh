#!/bin/bash

prefix=$1
datadir=$2
outdir=$3
chain=$4

mkdir -p $outdir

qsub -N ${chain}_$prefix processCDR3s_job.sh $prefix $datadir $outdir $chain
