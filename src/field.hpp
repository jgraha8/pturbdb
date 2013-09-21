#include "derivative.hpp"
#include "mpi_topology.hpp"

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

	int dims_[FIELD_NDIMS]; // Dimensions of global domain
	int dims_local_[FIELD_NDIMS]; // Dimensions of local domain
	int dims_operation_[FIELD_NDIMS]; // Dimension of "operations" domain

	int offset_local_[FIELD_NDIMS]; // Offset of local domain with respect to global domain
	int offset_operation_[FIELD_NDIMS]; // Offset of operation domain with respect to local domain

	int periodic_[FIELD_NDIMS]; // Settings for whether the domain in each field
	// direction is periodic or not
	int operator_order_; // Order of derivative operations
	int rind_size_; // Size of the MPI overlap rind

	// Decomposition
	FieldDecomp_t field_decomp_;

	// MPI Topology struct
	MpiTopology_t *mpi_topology_;

	// FiniteDiff Class 
	FiniteDiff *finite_diff_;
 

	// Grid pointers; these are not allocated but must be set with setGrid
	double *x_local_, *y_local_, *z_local_;

	bool synchronized_;

public:
	double *data_local; // Data of local domain

public:

	// Constructor
	Field(){};
	Field( const int *dims, FieldDecomp_t field_decomp, const int *periodic, int operator_order );
	Field(Field &g, bool copy_data_local);
	Field( Field &g );
	// Deconstructor
	~Field();

	// Derivative functions
	void finiteDiffInit();
    
	// Get the size of the field
	long           getSize();
	long           getSizeLocal();
	long           getSizeOperation();
	long           getSizeRind( int dim, int location );
	MpiTopology_t *getMpiTopology();
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

	long index( int i, int j, int k );
	long indexLocal( int i, int j, int k );
	long indexOperation( int i, int j, int k );
	long indexOperationToLocal( int i, int j, int k );

	void setGridLocal( double *x_local, double *y_local, double *z_local );
	void setDataOperation( const float *data_operation);
	void setDataOperation( const double *data_operation);

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
    
	void addDataOperation( const float *data_operation);
	void addDataOperation( const double *data_operation);
	void mulDataOperation( double scalar);

   
protected:

	void FieldInit( const int *dims, FieldDecomp_t field_decomp, const int *periodic, int operator_order );
	void FieldCopy( Field &g, bool copy_data_local );

	int *computeMpiTopologyDims( int nproc, int mpi_decomp_ndims );

	void assignMpiTopology();
	void assignDimsAndOffsets();

	bool hasRind(int dim, int location);
	double *createRindBuffer( int dim, int location );
	void packRindBuffer( int dim, int location, double *rind_buffer );
	void unpackRindBuffer( int dim, int location, double *rind_buffer );
	void synchronize();
	void synchronizeDimension(int coord_dim);

	void dndxn( void (FiniteDiff::*dd)(int, int, double *, int, double *), Field &a );
	void dndyn( void (FiniteDiff::*dd)(int, int, double *, int, double *), Field &a );
	void dndzn( void (FiniteDiff::*dd)(int, int, double *, int, double *), Field &a );
   
	// Getters/setters
	void setSynchronized( bool synchronized ); // Making protected since it should only be set by base class or sub-classes.

}; 

}
#endif
