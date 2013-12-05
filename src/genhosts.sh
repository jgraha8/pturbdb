#!/bin/bash

#NODES=( "dsp033" "dsp034" "dsp035" "dsp040" )
#NODES=( "dsp034" "dsp035" "dsp040" )
NODES=( "dsp034" )
PPN="$1"

echo -n > hostfile
for n in ${NODES[@]}
do 
    for ((p=0; p<$PPN; p++)) 
    do
	echo $n >> hostfile
    done
done
