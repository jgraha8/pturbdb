#include "Fields.hpp"
#include "mpi.h"

#ifndef MPIFIELDS_H
#define MPIFIELDS_H

namespace MPIFields 
{
  
  typedef enum {
    SLAB_DECOMP,
    PENCIL_DECOMP,
    CUBE_DECOMP
  } MPIDecomp_t;

  class MPIField: public Fields::Field // Inherits Field class
  {

  public: 

    MPI_Comm comm;
    int nproc, rank;
    MPIDecomp_t decomp;

    int *dims_global; // Global dimension

  public:
  
    // Constructor
    MPIField(){};
    MPIField( MPI_Comm _comm, int _nproc, int _rank, MPIDecomp_t _decomp, int _ndim, const int *_dims_global);
    MPIField( const MPIField &g );
    // Deconstructor
    ~MPIField();

    long getFieldAddr();
    long getSizeGlobal();

  private:
    
    void MPIFieldInit(MPI_Comm _comm, int _nproc, int _rank, MPIDecomp_t _decomp, int _ndim, const int *_dims);
    int *domainDecomp( int _nproc, int _rank, MPIDecomp_t _decomp, int _ndim, const int *_dims );

  }; 

}
#endif
