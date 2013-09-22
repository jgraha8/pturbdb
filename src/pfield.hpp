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
	PField &operator=( const PField &a );
	PField &operator=( double c );
	PField &operator+=( const PField &a );
	PField &operator+=( double c );
	PField &operator-=( const PField &a );
	PField &operator-=( double c );
	PField &operator*=( const PField &a );
	PField &operator*=( double c ) ;
	PField &operator/=( const PField &a );
	PField &operator/=( double c );

	PField &add( PField &a, PField &b );
	PField &sub( PField &a, PField &b );
	PField &mul( PField &a, PField &b );
	PField &div( PField &a, PField &b );

	void ddx( PField &a );
	void ddy( PField &a );
	void ddz( PField &a );

	void d2dx2( PField &a );
	void d2dy2( PField &a );
	void d2dz2( PField &a );
	void d2dxy( PField &a );
	void d2dxz( PField &a );
	void d2dyz( PField &a );
    
	void addDataOperation( const float *data_operation);
	void addDataOperation( const double *data_operation);
	void mulDataOperation( double scalar);

   
protected:

	void PFieldInit( const int *dims, FieldDecomp_t field_decomp, const int *periodic, int operator_order );
	void PFieldCopy( PField &g, bool copy_data_local );

	int *computeMpiTopologyDims( int nproc, int mpi_decomp_ndims );

	void assignMpiTopology();
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
