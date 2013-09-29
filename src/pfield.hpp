#include "derivative.hpp"
#include "mpi_topology.hpp"

#define PFIELD_NDIMS 3

#ifndef PFIELD_H
#define PFIELD_H

namespace pturbdb {

typedef enum {
	FIELD_DECOMP_SLAB,
	FIELD_DECOMP_PENCIL,
	FIELD_DECOMP_CUBE
} FieldDecomp_t;
  
class PField {

private:

	int dims_[PFIELD_NDIMS]; // Dimensions of global domain
	int dims_local_[PFIELD_NDIMS]; // Dimensions of local domain
	int dims_operation_[PFIELD_NDIMS]; // Dimension of "operations" domain

	int offset_local_[PFIELD_NDIMS]; // Offset of local domain with respect to global domain
	int offset_operation_[PFIELD_NDIMS]; // Offset of operation domain with respect to local domain

	int periodic_[PFIELD_NDIMS]; // Settings for whether the domain in each field
	// direction is periodic or not
	int operator_order_; // Order of derivative operations
	int rind_size_; // Size of the MPI overlap rind

	// Decomposition
	FieldDecomp_t field_decomp_;

	// MPI Topology struct
	MPITopology_t *mpi_topology_;

	// FiniteDiff Class 
	FiniteDiff *finite_diff_;
 

	// Grid pointers; these are not allocated but must be set with setGrid
	double *x_local_, *y_local_, *z_local_;

	bool synchronized_;

public:
	double *data_local; // Data of local domain

public:

	// Constructor
	PField(){};
	PField( const int *dims, FieldDecomp_t field_decomp, const int *periodic, int operator_order );
	PField(PField &g, bool copy_data_local);
	PField( PField &g );
	// Deconstructor
	~PField();

	// Derivative functions
	void finiteDiffInit();
    
	// Get the size of the field
	long           getSize();
	long           getSizeLocal();
	long           getSizeOperation();
	long           getSizeRind( int dim, int location );
	MPITopology_t *getMPITopology();
	FiniteDiff    *getFiniteDiff();
	FieldDecomp_t  getFieldDecomp();
	int           *getFieldPeriodic();
	int            getOperatorOrder();
	int           *getDims();
	int           *getDimsLocal();
	int           *getDimsOperation();
	int           *getOffsetLocal();
	int           *getOffsetOperation();
	int            getRindSize();
	bool           getSynchronized();

	double *getXLocal();
	double *getYLocal();
	double *getZLocal();	

	// Data getters; these require a buffer of size given by getDimsOperation()
	void getXOperation( double *x );
	void getYOperation( double *y );
	void getZOperation( double *z );

	// Data getters; these require a buffer of size given by getSizeOperation()
	void getDataOperation( float *a );
	void getDataOperation( double *a );

	long index( int i, int j, int k );
	long indexLocal( int i, int j, int k );
	long indexOperation( int i, int j, int k );
	long indexOperationToLocal( int i, int j, int k );

	void setGridLocal( double *x_local, double *y_local, double *z_local );
	void setDataOperation( const float *a ); // Perform same task as assignment operator
	void setDataOperation( const double *a ); // Performs same task as assignment operator

	// Assignment operator
	PField &operator=( const PField &a );
	PField &operator=( const double *a ); // Performs same task as setDataOperation
	PField &operator=( const float *a );  // Performs same task as setDataOperation
	PField &operator=( double c );

	PField &operator+=( const PField &a );
	PField &operator+=( const double *a );
	PField &operator+=( const float *a );
	PField &operator+=( double c );

	PField &operator-=( const PField &a );
	PField &operator-=( const double *a );
	PField &operator-=( const float *a );
	PField &operator-=( double c );

	PField &operator*=( const PField &a );
	PField &operator*=( const double *a );
	PField &operator*=( const float *a );
	PField &operator*=( double c ) ;

	PField &operator/=( const PField &a );
	PField &operator/=( const double *a );
	PField &operator/=( const float *a );
	PField &operator/=( double c );

	PField &add( PField &a, PField &b );
	PField &sub( PField &a, PField &b );
	PField &mul( PField &a, PField &b );
	PField &div( PField &a, PField &b );

	PField &ddx( PField &a );
	PField &ddy( PField &a );
	PField &ddz( PField &a );

	PField &d2dx2( PField &a );
	PField &d2dy2( PField &a );
	PField &d2dz2( PField &a );
	PField &d2dxy( PField &a );
	PField &d2dxz( PField &a );
	PField &d2dyz( PField &a );
    
protected:

	virtual void PFieldInit( const int *dims, FieldDecomp_t field_decomp, const int *periodic, int operator_order );
	void PFieldCopy( PField &g, bool copy_data_local );

	int *computeMPITopologyDims( int nproc, int mpi_decomp_ndims );

	void assignMPITopology();
	void assignDimsAndOffsets();

	bool hasRind(int dim, int location);
	double *createRindBuffer( int dim, int location );
	void packRindBuffer( int dim, int location, double *rind_buffer );
	void unpackRindBuffer( int dim, int location, double *rind_buffer );
	void synchronize();
	void synchronizeDimension(int coord_dim);

	void dndxn( void (FiniteDiff::*dd)(int, int, double *, int, double *), PField &a );
	void dndyn( void (FiniteDiff::*dd)(int, int, double *, int, double *), PField &a );
	void dndzn( void (FiniteDiff::*dd)(int, int, double *, int, double *), PField &a );
   
	// Getters/setters
	void setSynchronized( bool synchronized ); // Making protected since it should only be set by base class or sub-classes.

}; 

}
#endif
