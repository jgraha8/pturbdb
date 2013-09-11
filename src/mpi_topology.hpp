#include "mpi.h"

#define MPI_TOPOLOGY_CART_REORDER 1

#ifndef MPI_TOPOLOGY_HPP
#define MPI_TOPOLOGY_HPP

typedef struct {
  MPI_Comm comm;
  int nproc, rank;
  int ndims;
  int *dims;
  int *periodic;
  int *coords;
  int *neighbor_next;
  int *neighbor_prev;
} MpiTopology_t;

#endif

MpiTopology_t *MpiTopologyNew( MPI_Comm comm, int ndims, int *dims, int *periodic );
void MpiTopologyDelete( MpiTopology_t **mpi_topology );
