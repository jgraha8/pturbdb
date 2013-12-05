#ifndef PFIELD_H
#define PFIELD_H


#include <algorithm>
#include <vector>
#include "derivative.hpp"
#include "mpi_topology.hpp"

#define PFIELD_NDIMS 3

#define PFIELD_LOOP_OPERATION(pfield) \
	for( int _i=0; _i<pfield->getDimsOperation()[0]; _i++ ) {	\
	        for( int _j=0; _j<pfield->getDimsOperation()[1]; _j++ ) { \
		        for( int _k=0; _k<pfield->getDimsOperation()[2]; _k++ ) { \
			        const size_t _index = pfield->indexOperation(_i,_j,_k); 

#define PFIELD_LOOP_OPERATION_TO_LOCAL(pfield) \
	for( int _i=0; _i<pfield->getDimsOperation()[0]; _i++ ) { \
	        for( int _j=0; _j<pfield->getDimsOperation()[1]; _j++ ) { \
		        for( int _k=0; _k<pfield->getDimsOperation()[2]; _k++ ) { \
			        const size_t _index = pfield->indexOperationToLocal(_i,_j,_k); 
#define PFIELD_LOOP_LOCAL(pfield) \
	for( int _i=0; _i<pfield->getDimsLocal()[0]; _i++ ) {	\
	        for( int _j=0; _j<pfield->getDimsLocal()[1]; _j++ ) { \
		        for( int _k=0; _k<pfield->getDimsLocal()[2]; _k++ ) { \
			        const size_t _index = pfield->indexLocal(_i,_j,_k); 

#define PFIELD_LOOP_END }}}

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
	virtual ~PField();

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
	const double  *getDataLocal()         const { return this->data_local;       }

	// Data getters; these require a buffer of size given by getDimsOperation()
	void getXOperation( double *x ) const;
	void getYOperation( double *y ) const;
	void getZOperation( double *z ) const;

	// Data getters; these require a buffer of size given by getSizeOperation()
	void getDataOperation( float *a ) const;
	void getDataOperation( double *a ) const;

	size_t index( const int &i, const int &j, const int &k ) const;
	size_t indexLocal( const int &i, const int &j, const int &k ) const;
	size_t indexOperation( const int &i, const int &j, const int &k ) const;
	size_t indexOperationToLocal( const int &i, const int &j, const int &k ) const;

	std::vector<int> ijk( const size_t &index ) const;
	std::vector<int> ijkLocal( const size_t &index ) const;
	std::vector<int> ijkOperation( const size_t &index ) const;

	bool inDomain( const int &i, const int &j, const int &k ) const;
	bool inDomainLocal( const int &i, const int &j, const int &k ) const;
	bool inDomainOperation( const int &i, const int &j, const int &k ) const;

	void setGridLocal( const double *x_local, const double *y_local, const double *z_local );
	void setDataLocal( const size_t &index, const double &a );
	void setDataOperation( const float *a ); // Perform same task as assignment operator
	void setDataOperation( const double *a ); // Performs same task as assignment operator
	void setDataOperation( const int &i, const int &j, const int &k, const double &a );

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

	PField &filter( const int &filter_width );

	void synchronize();
    
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
	void synchronizeDimension(int coord_dim);

	void dndxn( void (FiniteDiff::*dd)(int, int, const double *, int, double *), PField &a );
	void dndyn( void (FiniteDiff::*dd)(int, int, const double *, int, double *), PField &a );
	void dndzn( void (FiniteDiff::*dd)(int, int, const double *, int, double *), PField &a );

	// The index i is for the operation domain
	inline const int filter_width_left_boundary( const int &dim, 
							      const int &i, 
							      const int &filter_width_half ) {
		const int i_global = i + this->offset_operation_[dim] + this->offset_local_[dim];	
		return std::min(i_global, filter_width_half );	
	}
	inline const int filter_width_right_boundary( const int &dim, 
							       const int &i, 
							       const int &filter_width_half ) {
		const int i_global = i + this->offset_operation_[dim] + this->offset_local_[dim];	
		return std::min(this->dims_[dim] - 1 - i_global, filter_width_half );	
	}
	inline const int filter_width_left_rind( const int &dim, 
							  const int &i, 
							  const int &filter_width_half ) {
		return filter_width_half;
	}
	inline const int filter_width_right_rind( const int &dim, 
							   const int &i, 
							   const int &filter_width_half ) {
		return filter_width_half;
	}

   
	// Getters/setters
	void setSynchronized( bool synchronized ); // Making protected since it should only be set by base class or sub-classes.

}; 

}
#endif
