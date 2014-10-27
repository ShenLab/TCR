#!/bin/bash
#$ -cwd
#$ -o log
#$ -e log

time1=$( date "+%s" )

echo [BEGIN] `date`
echo [MACHINE] `hostname`

BWA=/ifs/scratch/c2b2/ys_lab/aps2157/software/bwa
SAMTOOLS=/ifs/scratch/c2b2/ys_lab/bg2178/shared/software/samtools-1.0/samtools

fqfile=$1
ref=$2
outdir=$3

echo $ref
echo $fqfile

name=`basename $fqfile`
$BWA bwasw -t 4 -z 3 -c 3 -r 1 -q 4 -w 100 $ref $fqfile > $outdir/$name.sam
$SAMTOOLS view -bS $outdir/$name.sam > $outdir/$name.bam
rm $outdir/$name.sam

time2=$( date "+%s" )
echo [deltat] $(( $time2 - $time1 ))

echo [END] `date`
