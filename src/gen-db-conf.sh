#!/bin/bash

DB_ROOT="/datascope/tdbchannel"

echo "${DB_ROOT}/grid.h5" > turbdb.conf

dt=0.0013
N_step_db=5;

index=0;
for ((i=1; i<=10; i++))
do
    for f in ${DB_ROOT}/tdb-chunk-$i/fields*.h5
    do 
	db_time=$(echo "scale=9; $((index-1))*${dt}*${N_step_db}" | bc)
	echo $db_time $f >> turbdb.conf
	((index++))
    done
done
