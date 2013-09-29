#include "mpi.h"

#define MPI_TOPOLOGY_CART_REORDER 1

#ifndef MPI_TOPOLOGY_HPP
#define MPI_TOPOLOGY_HPP

typedef struct {
  MPI_Comm comm;
  int nproc, rank;
  int ndims;
  int *dims;
  int *coords;
  int *neighbor_next;
  int *neighbor_prev;
} MPITopology_t;

#endif

MPITopology_t *MPITopologyNew( MPI_Comm comm, int ndims, int *dims, int *periodic );
void MPITopologyDelete( MPITopology_t **mpi_topology );
