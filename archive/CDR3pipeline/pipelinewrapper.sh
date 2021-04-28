#!/bin/bash

prefix=$1
datadir=$2
outdir=$3
chain=$4

if [ $chain == A ]; then
	code=doTRA.sh
elif [ $chain == B ]; then
	code=doTRB.sh
else
	echo "Improper input chain = $chain"
fi

for i in $datadir/$prefix*.bam;
do 	
	qsub -N ${chain}_`basename $i` $code $i $outdir 
done


