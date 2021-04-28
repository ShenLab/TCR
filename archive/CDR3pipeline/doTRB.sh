#!/bin/bash
#$ -cwd
#$ -l mem=10G,time=48::
#$ -o log
#$ -e log

time1=$( date "+%s" )

echo [BEGIN] `date`
echo [MACHINE] `hostname`
                              
sample=$1 # sample name (bam file)
TRBoutdir=$2/`basename $1`.Bchain # output directory

SAMTOOLS=/ifs/scratch/c2b2/ys_lab/bg2178/shared/software/samtools-1.0/samtools
mkdir -p $TRBoutdir

$SAMTOOLS view -F 4 $sample > `basename $sample .bam`.sam
perl ReadSam.CDR3.TRB.pl $TRBoutdir `basename $sample .bam`.sam
perl ReadCDR3.VJ.TRB.pl $TRBoutdir `basename $a .bam`.sam.CDR3.VJ.seq
perl CategorizeCDR3.TRB.pl $TRBoutdir `basename $sample .bam`.sam.CDR3.VJ.seq.prot.txt
rm `basename $sample .bam`.sam


time2=$( date "+%s" )
echo [deltat] $(( $time2 - $time1 ))

echo [END] `date`
