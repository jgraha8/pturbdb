#!/bin/bash

EXEC=channel_q_field
#SRCS=( derivative.cpp pfield.cpp pfield_math.cpp pturbdb_field.cpp mpi_topology.cpp turbdb_analysis.cpp )
SRCS=( derivative.cpp pfield.cpp pfield_math.cpp pturbdb_field.cpp mpi_topology.cpp vortex_tracking.cpp ${EXEC}.cpp )
# Generate the object list
OBJS=$(for s in ${SRCS[@]}; do echo ${s%.*}.o; done)

INCPATH="-I../autofd -I${ESIO_ROOT}/include -I${LAPACK_ROOT}/include"
LIBPATH="-L../autofd -L../matsolv -L${ESIO_ROOT}/lib -L${HDF5_ROOT}/lib -L${LAPACK_ROOT}/lib"
LIBS="-lautofd -lmatsolv -lgfortran -lesio -lhdf5_hl -lhdf5 -lz -llapacke -llapack"

CFLAGS="-O2 -Wall -g -DBOUNDS_CHECK"

for s in ${SRCS[@]}
do 
    o=${s%.*}.o
    rm -f $o
    echo mpicxx $CFLAGS $INCPATH -c -o $o $s
    mpicxx $CFLAGS $INCPATH -c -o $o $s
done

# Link everything
#echo "mpicxx $LIBPATH -o speed3 $OBJS $LIBS"
#mpicxx $LIBPATH -o speed3 $OBJS $LIBS

echo "mpicxx $LIBPATH -o $EXEC $OBJS $LIBS"
mpicxx $LIBPATH -o $EXEC $OBJS $LIBS
