#!/bin/bash

EXEC=speed3
SRCS=( derivative.cpp field.cpp mpi_topology.cpp ) #speed3.cpp )
# Generate the object list
OBJS=$(for s in ${SRCS[@]}; do echo ${s%.*}.o; done)

INCPATH="-I../../autofd"
LIBPATH="-L../../autofd -L../../matsolv"
LIBS="-lautofd -lmatsolv -lgfortran"

CFLAGS="-O0 -Wall -g -DBOUNDS_CHECK"

for s in ${SRCS[@]}
do
    o=${s%.*}.o
    echo mpicxx $CFLAGS $INCPATH -c -o $o $s
    mpicxx $CFLAGS $INCPATH -c -o $o $s
done

# # Link everything
# echo "mpicxx $LIBPATH -o speed3 $OBJS $LIBS"
# mpicxx $LIBPATH -o speed3 $OBJS $LIBS
