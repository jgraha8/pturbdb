#!/bin/bash

EXEC=speed3
SRCS=( derivative.cpp field.cpp turbdb_field.cpp mpi_topology.cpp turbdb_analysis.cpp )
# Generate the object list
OBJS=$(for s in ${SRCS[@]}; do echo ${s%.*}.o; done)

INCPATH="-I../autofd -I${ESIO_ROOT}/include"
LIBPATH="-L../autofd -L../matsolv -L${ESIO_ROOT}/lib -L${HDF5_ROOT}/lib"
LIBS="-lautofd -lmatsolv -lgfortran -lesio -lhdf5_hl -lhdf5 -lz"

CFLAGS="-O0 -Wall -g -DBOUNDS_CHECK"

for s in ${SRCS[@]}
do
    o=${s%.*}.o
    echo mpicxx $CFLAGS $INCPATH -c -o $o $s
    mpicxx $CFLAGS $INCPATH -c -o $o $s
done

# Link everything
#echo "mpicxx $LIBPATH -o speed3 $OBJS $LIBS"
#mpicxx $LIBPATH -o speed3 $OBJS $LIBS

echo "mpicxx $LIBPATH -o turbdb_analysis $OBJS $LIBS"
mpicxx $LIBPATH -o turbdb_analysis $OBJS $LIBS
