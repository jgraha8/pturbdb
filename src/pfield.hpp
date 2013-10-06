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
	const double *x_local_, *y_local_, *z_local_;

	bool synchronized_;

public:
	double *data_local; // Data of local domain

public:

	// Constructor
	PField(){};
	PField( const int *dims, FieldDecomp_t field_decomp, const int *periodic, int operator_order );
	PField( const PField &g, bool copy_data_local);
	PField( const PField &g );
	// Deconstructor
	~PField();

	// Derivative functions
	void finiteDiffInit();
    
	// Get the size of the field
	size_t getSize() const;
	size_t getSizeLocal() const;
	size_t getSizeOperation() const;
	size_t getSizeRind( int dim, int location ) const;


	FieldDecomp_t  getFieldDecomp()       const { return this->field_decomp_;     };
	const MPITopology_t *getMPITopology() const { return this->mpi_topology_;     };
	const FiniteDiff    *getFiniteDiff()  const { return this->finite_diff_;      };
	const double        *getXLocal()      const { return this->x_local_;          };
	const double        *getYLocal()      const { return this->y_local_;          };
	const double        *getZLocal()      const { return this->z_local_;          };
	const int     *getFieldPeriodic()     const { return this->periodic_;         };
	int            getOperatorOrder()     const { return this->operator_order_;   };
	const int     *getDims()              const { return this->dims_;             };
	const int     *getDimsLocal()         const { return this->dims_local_;       };    
	const int     *getDimsOperation()     const { return this->dims_operation_;   };
	const int     *getOffsetLocal()       const { return this->offset_local_;     };
	const int     *getOffsetOperation()   const { return this->offset_operation_; };
	int            getRindSize()          const { return this->rind_size_;        };
	bool           getSynchronized()      const { return this->synchronized_;     };

	// Data getters; these require a buffer of size given by getDimsOperation()
	void getXOperation( double *x ) const;
	void getYOperation( double *y ) const;
	void getZOperation( double *z ) const;

	// Data getters; these require a buffer of size given by getSizeOperation()
	void getDataOperation( float *a ) const;
	void getDataOperation( double *a ) const;

	size_t index( int i, int j, int k ) const;
	size_t indexLocal( int i, int j, int k ) const;
	size_t indexOperation( int i, int j, int k ) const;
	size_t indexOperationToLocal( int i, int j, int k ) const;

	void setGridLocal( const double *x_local, const double *y_local, const double *z_local );
	void setDataOperation( const float *a ); // Perform same task as assignment operator
	void setDataOperation( const double *a ); // Performs same task as assignment operator

	// Assignment operator
	PField &operator=( float c);
	PField &operator=( double c );
	PField &operator=( const float *a );  // Performs same task as setDataOperation
	PField &operator=( const double *a ); // Performs same task as setDataOperation
	PField &operator=( const PField &a );

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
	PField &operator*=( float c ) ;
	PField &operator*=( double c ) ;

	PField &operator/=( const PField &a );
	PField &operator/=( const double *a );
	PField &operator/=( const float *a );
	PField &operator/=( double c );

	PField &add( const PField &a, const PField &b );
	PField &sub( const PField &a, const PField &b );
	PField &mul( const PField &a, const PField &b );
	PField &div( const PField &a, const PField &b );

	PField &sqrt( const PField &a );

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
	void PFieldCopy( const PField &g, bool copy_data_local );

	int *computeMPITopologyDims( int nproc, int mpi_decomp_ndims ) const;

	void assignMPITopology();
	void assignDimsAndOffsets();

	bool hasRind(int dim, int location) const;
	double *createRindBuffer( int dim, int location ) const;
	void packRindBuffer( int dim, int location, double *rind_buffer ) const;
	void unpackRindBuffer( int dim, int location, const double *rind_buffer );
	void synchronize();
	void synchronizeDimension(int coord_dim);

	void dndxn( void (FiniteDiff::*dd)(int, int, const double *, int, double *), PField &a );
	void dndyn( void (FiniteDiff::*dd)(int, int, const double *, int, double *), PField &a );
	void dndzn( void (FiniteDiff::*dd)(int, int, const double *, int, double *), PField &a );
   
	// Getters/setters
	void setSynchronized( bool synchronized ); // Making protected since it should only be set by base class or sub-classes.

}; 

}
#endif
