#include "derivative.hpp"
#incllude "mpi_topology.hpp"

#define FIELD_NDIMS 3

#ifndef FIELD_H
#define FIELD_H

namespace pturb_fields {

  typedef enum {
    FIELD_DECOMP_SLAB,
    FIELD_DECOMP_PENCIL,
    FIELD_DECOMP_CUBE
  } FieldDecomp_t;
  
  class Field {

  private:

    int *dims_; // Dimensions of global domain
    int *dims_local_; // Dimensions of local domain
    int *dims_operation_; // Dimension of "operations" domain

    int *offset_local_; // Offset of local domain with respect to global domain
    int *offset_operation_; // Offset of operation domain with respect to local domain

    // Grid pointers; these are not allocated but must be set with setGrid
    double *x_local_, *y_local_, *z_local_;

    // FiniteDiff Class 
    FiniteDiff *finite_diff_;
 
    // Decomposition
    FieldDecomp_t field_decomp_;

    // MPI Topology struct
    MpiTopology_t *mpi_topology_;

  public:
    double *data_local; // Data of local domain

  public:

    // Constructor
    Field(){};
    Field( int ndim, int *dims, int *periodic, FieldDecomp_t field_decomp, int operator_order );
    Field( const Field &g );
    // Deconstructor
    ~Field();


    // Derivative functions
    void derivFDInit( int order );
    
    // Get the size of the field
    long getSize();
    MPI_Comm getMpiComm();
    FieldDecomp_t getFieldDecomp();
    int *getDims();

    long index( int i, int j, int k );
    long indexLocal( int i, int j, int k );
    long indexOperation( int i, int j, int k );
    long indexOperationToLocal( int i, int j, int k );

    void assignGridLocal( double *x_local, double *y_local, double *z_local );

    // Assignment operator
    Field &operator=( const Field &a );
    Field &operator+=( const Field &a );
    Field &operator-=( const Field &a );
    Field &operator*=( const Field &a );
    Field &operator/=( const Field &a );

    void add( Field &a, Field &b );
    void sub( Field &a, Field &b );
    void mul( Field &a, Field &b );
    void div( Field &a, Field &b );

    void ddx( Field &a );
    void ddy( Field &a );
    void ddz( Field &a );

    void d2dx2( Field &a );
    void d2dy2( Field &a );
    void d2dz2( Field &a );
    void d2dxy( Field &a );
    void d2dxz( Field &a );
    void d2dyz( Field &a );
    
   
  private:

    void FieldInit( int *dims );
    void FieldInit( int ndim, int *dims );

    int *domainDecomp( int _nproc, int _rank, MPIDecomp_t _decomp, int _ndim, const int *_dims );

    void dndxn( void (FiniteDiff::*dd)( int, double *, double *), Field &a );
    void dndyn( void (FiniteDiff::*dd)( int, double *, double *), Field &a );
    void dndzn( void (FiniteDiff::*dd)( int, double *, double *), Field &a );
   

  }; 

}

#endif
