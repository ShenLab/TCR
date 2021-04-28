#!/bin/bash
#$ -cwd
#$ -l mem=5G,time=2::
#$ -o log
#$ -e log


prefix=$1
datadir=$2
outdir=$3
chain=$4


# get all CDR3 sequences
for i in $datadir/$prefix*${chain}chain
do
	cat $i/*.seq.prot.txt | cut -f4,5,12,15 | sort -k3 | uniq -c | awk '{print $2"."$3"\t"$1"\t"$2"\t"$3"\t"$4"\t"$5"\t"length($4)}' > $outdir/$prefix.${chain}chain.full 
done


# get productive CDR3 sequences
infile=$outdir/$prefix.${chain}chain.full
echo -e "VJCombine\tCopy\tVGeneName\tJGeneName\tCDR3\tSequence\tLengthOfCDR3" >$outdir/`basename $infile .full`.productive.tsv ; cat $infile |grep -v Out | grep -v "*" | awk '{if(NF==7) print $0}' | sort -k2 -nr >> $outdir/`basename $infile .full`.productive.tsv;


# final error correction step
infile=$outdir/`basename $infile .full`.productive.tsv
Jmotif=human${chain}.Jmotif

echo -e "VJCombine\tCopy\tVGeneName\tJGeneName\tCDR3\tSequence" > $outdir/`basename $infile .productive.tsv`.final.tsv
sed 1d $infile |  python fixNTs.py $Jmotif  |  awk -f merge.awk | tr @ \\t | sort -rnk6 | awk '{print $1"\t"$6"\t"$2"\t"$3"\t"$4"\t"$5}' >> $outdir/`basename $infile .productive.tsv`.final.tsv
