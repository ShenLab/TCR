#!/bin/bash
#$ -cwd
#$ -l mem=10G,time=48::
#$ -o log
#$ -e log

# Main pipeline to run the alpha chain analysis

time1=$( date "+%s" )

echo [BEGIN] `date`
echo [MACHINE] `hostname`

sample=$1 # sample name (bam file)
TRAoutdir=$2/`basename $1`.Achain # output directory


SAMTOOLS=/ifs/scratch/c2b2/ys_lab/bg2178/shared/software/samtools-1.0/samtools
mkdir -p $TRAoutdir

$SAMTOOLS view -F 4 $sample > `basename $sample .bam`.sam
perl ReadSam.CDR3.TRA.pl $TRAoutdir `basename $sample .bam`.sam # find VJ cassette and CDR3 region
perl ReadCDR3.VJ.TRA.pl $TRAoutdir `basename $sample .bam`.sam.CDR3.VJ.seq # Get reading frame and translate
perl CategorizeCDR3.TRA.pl $TRAoutdir `basename $sample .bam`.sam.CDR3.VJ.seq.prot.txt # collect output
rm `basename $sample .bam`.sam

time2=$( date "+%s" )
echo [deltat] $(( $time2 - $time1 ))

echo [END] `date`
