#include <cstring>
#include "mpi_topology.hpp"

MPITopology_t *MPITopologyNew( MPI_Comm comm, int ndims, int *dims, int *periodic )
{

  MPITopology_t *mpi_topology = new MPITopology_t;
  
  mpi_topology->ndims         = ndims;
  mpi_topology->dims          = new int[ndims];
  mpi_topology->coords        = new int[ndims];
  mpi_topology->neighbor_next = new int[ndims];
  mpi_topology->neighbor_prev = new int[ndims];

  memcpy(mpi_topology->dims,dims,sizeof(dims[0])*ndims);

  MPI_Cart_create( comm, 
		   mpi_topology->ndims, 
		   mpi_topology->dims, 
		   periodic,
		   MPI_TOPOLOGY_CART_REORDER, 
		   &mpi_topology->comm );


  MPI_Comm_rank( mpi_topology->comm, &mpi_topology->rank );
  MPI_Comm_size( mpi_topology->comm, &mpi_topology->nproc );

  MPI_Cart_coords( mpi_topology->comm, 
		   mpi_topology->rank, 
		   mpi_topology->ndims, 
		   mpi_topology->coords );
 
  // Get the process neighbors in each direction
  for( int n=0; n<mpi_topology->ndims; n++ ) {
    MPI_Cart_shift( mpi_topology->comm, n, 1, &mpi_topology->neighbor_prev[n], &mpi_topology->neighbor_next[n] );
  }

  
  return mpi_topology;

};

void MPITopologyDelete( MPITopology_t **mpi_topology ) 
{

  delete [] (*mpi_topology)->dims; 
  delete [] (*mpi_topology)->coords; 
  delete [] (*mpi_topology)->neighbor_next; 
  delete [] (*mpi_topology)->neighbor_prev; 
  delete *mpi_topology;

  *mpi_topology = NULL;

}
