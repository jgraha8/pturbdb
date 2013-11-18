#include <iostream>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <cstdio>
#include "pfield.hpp"

using namespace std;

namespace pturbdb {

////////////////////////////////////////////////////////////////////////////////
/// CONSTRUCTORS
////////////////////////////////////////////////////////////////////////////////

// PField class member functions
PField::PField(const int *dims, FieldDecomp_t field_decomp, const int *periodic,
	     int operator_order)
{
	PFieldInit(dims, field_decomp, periodic, operator_order);

#ifdef VERBOSE
	printf("%d: PField::PField: constructed field\n", this->getMPITopology()->rank);
#endif

}

PField::PField(const PField &g, bool copy_data_local)
{

	this->PFieldCopy( g, copy_data_local );

#ifdef VERBOSE
	printf("%d: PField::PField: constructed field from copy\n", this->getMPITopology()->rank);
#endif

}

PField::PField(const PField &g)
{

	this->PFieldCopy( g, true );

#ifdef VERBOSE
	printf("%d: PField::PField: constructed field from copy\n", this->getMPITopology()->rank);
#endif

}

////////////////////////////////////////////////////////////////////////////////
/// DECONSTRUCTORS
////////////////////////////////////////////////////////////////////////////////

PField::~PField(void)
{

#ifdef VERBOSE
	printf("%d: PField::~PField: deconstructed field\n", this->mpi_topology_->rank);
#endif

	if( this->mpi_topology_ != NULL ) MPITopologyDelete( &this->mpi_topology_ );
	if( this->finite_diff_ != NULL ) delete this->finite_diff_;
       	delete[] this->data_local;



}

////////////////////////////////////////////////////////////////////////////////
/// INITIALIZERS
////////////////////////////////////////////////////////////////////////////////

// PField class initializer
// Sets:
//   ndims
// Initilizes:
//   dims
//   data
/********************************************************************/
void PField::PFieldInit(const int *dims, FieldDecomp_t field_decomp,
		      const int *periodic, int operator_order)
/********************************************************************/
{

	// Set the global dimensions
	memcpy(this->dims_, dims, sizeof(*this->dims_) * PFIELD_NDIMS);
	this->field_decomp_ = field_decomp;
	memcpy(this->periodic_, periodic, sizeof(*this->periodic_) * PFIELD_NDIMS);
	this->operator_order_ = operator_order;
	// Rind size to support finite differencing at interprocessor boundaries
	this->rind_size_ = operator_order / 2 + operator_order % 2; // Second part adds 1 if odd operator order

	// Compute the MPI topology for the given field decomposition
	this->assignMPITopology();

	// Compute the local and operation dims and offsets
	this->assignDimsAndOffsets();

	// Now that we have set dims_local_ we may get the total number of elements
	this->data_local = new double[this->getSizeLocal()];

	// Initialize pointers to null
	this->x_local_ = NULL;
	this->y_local_ = NULL;
	this->z_local_ = NULL;
	this->finite_diff_ = NULL;

	// Initalize state variables
	this->synchronized_ = false;

#ifdef VERBOSE
	for( int n=0; n<this->mpi_topology_->nproc; n++) {
		if( n == this->mpi_topology_->rank ) {
			printf("%d: PField::PFieldInit\n", n);
			printf("    dims             = %d, %d\n", dims_[0], dims_[1]);
			printf("    dims_local       = %d, %d\n", dims_local_[0], dims_local_[1]);
			printf("    dims_operation   = %d, %d\n", dims_operation_[0], dims_operation_[1]);
			printf("    offset_local     = %d, %d\n", offset_local_[0], offset_local_[1]);
			printf("    offset_operation = %d, %d\n", offset_operation_[0], offset_operation_[1]);
			printf("    periodic         = %d, %d\n", periodic_[0], periodic_[1]);
			printf("    operator_order   = %d\n", operator_order_);
			printf("    rind_size        = %d\n", rind_size_);
		}
		MPI_Barrier(this->mpi_topology_->comm);
	}
#endif

}

/*
 * Worker for the PField copy constructor
 */
void PField::PFieldCopy( const PField &g, bool copy_data_local )
{

	// We will make a copy of the input field
	// Copy dimension data
	memcpy(this->dims_,             g.getDims(),            sizeof(this->dims_[0])*PFIELD_NDIMS);
	memcpy(this->dims_local_,       g.getDimsLocal(),       sizeof(this->dims_local_[0])*PFIELD_NDIMS);
	memcpy(this->dims_operation_,   g.getDimsOperation(),   sizeof(this->dims_operation_[0])*PFIELD_NDIMS);

	// Copy offsets
	memcpy(this->offset_local_,     g.getOffsetLocal(),     sizeof(this->offset_local_[0])*PFIELD_NDIMS);
	memcpy(this->offset_operation_, g.getOffsetOperation(), sizeof(this->offset_operation_[0])*PFIELD_NDIMS);
    
	// Copy periodic data
	memcpy(this->periodic_,         g.getFieldPeriodic(),   sizeof(this->periodic_[0])*PFIELD_NDIMS);

	// Copy scalars
	this->operator_order_ = g.getOperatorOrder();
	this->rind_size_      = g.getRindSize();
	this->field_decomp_   = g.getFieldDecomp();

	// Copy external pointers
	this->x_local_        = g.getXLocal();
	this->y_local_        = g.getYLocal();
	this->z_local_        = g.getZLocal();

	// Assign the MPI topology. Cannot copy the pointer value
	// since if g is deleted we lose the allocated MPI topology
	// struct.
	this->assignMPITopology();
	this->finite_diff_ = NULL;
	// Initialize finite difference class
	if( this->x_local_ != NULL && this->y_local_ != NULL && this->z_local_ != NULL ) {
		this->finiteDiffInit();
	} 

	// Copy synchronization state
	this->synchronized_   = g.getSynchronized();
 
	// Allocate and copy the data on the local domain
	this->data_local = new double[this->getSizeLocal()];

	if( copy_data_local )
		memcpy(this->data_local, g.data_local, sizeof(*this->data_local) * this->getSizeLocal());

}

//
// Initialize the finite difference derivatives class. Requires that
// setGridLocal be called from the application before this can be called.
//
/********************************************************************/
void PField::finiteDiffInit()
/********************************************************************/
{


#ifdef VERBOSE
	printf("%d: PField::finiteDiffInit: entering\n", this->mpi_topology_->rank);
#endif

	int ierr=0;

	if (this->finite_diff_ != NULL) {
		printf("%d: PField::finiteDiffInit: class instance fd already created--skipping.\n", this->mpi_topology_->rank);
		//MPI_Abort(this->mpi_topology_->comm, ierr);
	}

	if (this->x_local_ == NULL || this->y_local_ == NULL || this->z_local_ == NULL) {
		printf("%d: PField::finiteDiffInit: requires setGridLocal be called first\n", this->mpi_topology_->rank);
		MPI_Abort(this->mpi_topology_->comm, ierr);
	}

	// First initialize the fd class instance
	this->finite_diff_ = new FiniteDiff(this->dims_local_[0], this->x_local_,
					    this->dims_local_[1], this->y_local_, 
					    this->dims_local_[2], this->z_local_, 
					    this->operator_order_);

#ifdef VERBOSE
	printf("%d: PField::finiteDiffInit: exiting\n", this->mpi_topology_->rank);
#endif

}

//////////////////////////////////////////////////////////////////////
/// GETTERS/SETTERS
//////////////////////////////////////////////////////////////////////

/********************************************************************/
size_t PField::getSize() const
/********************************************************************/
{
	size_t N = 1;
	for (int n = 0; n < PFIELD_NDIMS; n++) {
		N *= this->dims_[n];
	}
	return N;
}

/********************************************************************/
size_t PField::getSizeLocal() const
/********************************************************************/
{
	size_t N = 1;
	for (int n = 0; n < PFIELD_NDIMS; n++) {
		N *= this->dims_local_[n];
	}
	return N;
}

/********************************************************************/
size_t PField::getSizeOperation() const
/********************************************************************/
{
	size_t N = 1;
	for (int n = 0; n < PFIELD_NDIMS; n++) {
		N *= this->dims_operation_[n];
	}
	return N;
}

/********************************************************************/
size_t PField::getSizeRind(int dim, int location) const
/********************************************************************/
{
	// If the dimension and location does not have a rind, return 0
	if ( !this->hasRind(dim, location) ) 
		return 0;

	if (dim == 0) {
		return (size_t) this->rind_size_ * (size_t) this->dims_operation_[1]
			* (size_t) this->dims_operation_[2];
	} else if (dim == 1) {
		return (size_t) this->dims_operation_[0] * (size_t) this->rind_size_
			* (size_t) this->dims_operation_[2];
	} else if (dim == 2) {
		return (size_t) this->dims_operation_[0] * (size_t) this->dims_operation_[1]
			* (size_t) this->rind_size_;
	} else {
		printf("%d: PField::getSizeRind: incorrect dimension number specified\n", 
		       this->mpi_topology_->rank);
		int ierr=0;
		MPI_Abort(this->mpi_topology_->comm, ierr);
		return 0; // This shouldn't be reached but put here to avoid compiler warning
	}
}


/*
 * Getters for data members
 */

/*
 * Gets the x data on the operation grid by copying to a provided data
 * buffer. The size of "x" must be of size getDimsOperation().
 */ 
void PField::getXOperation( double *x) const
{
	for (int i = 0; i < this->dims_operation_[0]; i++)
		x[i] = this->x_local_[i + this->offset_operation_[0]];
}
/*
 * Gets the y data on the operation grid by copying to a provided data
 * buffer. The size of "y" must be of size getDimsOperation().
 */ 
void PField::getYOperation( double *y) const
{
	for (int j = 0; j < this->dims_operation_[1]; j++)
		y[j] = this->y_local_[j + this->offset_operation_[1]];
}
/*
 * Gets the z data on the operation grid by copying to a provided data
 * buffer. The size of "z" must be of size getDimsOperation().
 */ 
void PField::getZOperation( double *z) const
{
	for (int k = 0; k < this->dims_operation_[2]; k++)
		z[k] = this->z_local_[k + this->offset_operation_[2]];
}


/*
 * Gets the data on the operation grid by copying to a provided data
 * buffer. The size of "a" must be of size getSizeOperation().
 */ 
void PField::getDataOperation( float *a) const
{
	size_t index_operation;

	index_operation=0;
	PFIELD_LOOP_OPERATION_TO_LOCAL(this)
		a[index_operation++] = (float)this->data_local[_index];
	PFIELD_LOOP_END
}

/*
 * Gets the data on the operation grid by copying to a provided data
 * buffer. The size of "a" must be of size getSizeOperation().
 */ 
void PField::getDataOperation( double *a) const
{
	size_t index_operation;

	index_operation=0;
	PFIELD_LOOP_OPERATION_TO_LOCAL(this)
		a[index_operation++] = this->data_local[_index];
	PFIELD_LOOP_END
}


// Array index for 3D indexing
/********************************************************************/
size_t PField::index( const int &i, const int &j, const int &k ) const
/********************************************************************/
{
	return ((size_t) this->dims_[1] * (size_t) i + (size_t) j) * (size_t) this->dims_[2]
		+ (size_t) k;
}

// Array index for 3D indexing
/********************************************************************/
size_t PField::indexLocal(const int &i, const int &j, const int &k) const
/********************************************************************/
{
	return ((size_t) this->dims_local_[1] * (size_t) i + (size_t) j)
		* (size_t) this->dims_local_[2] + (size_t) k;
}

// Array index for 3D indexing
/********************************************************************/
size_t PField::indexOperation(const int &i, const int &j, const int &k) const
/********************************************************************/
{
	return ((size_t) this->dims_operation_[1] * (size_t) i + (size_t) j)
		* (size_t) this->dims_operation_[2] + (size_t) k;
}

// Array index for 3D indexing
/********************************************************************/
size_t PField::indexOperationToLocal(const int &i, const int &j, const int &k) const
/********************************************************************/
{
	const size_t ii = (size_t) i + (size_t) this->offset_operation_[0];
	const size_t jj = (size_t) j + (size_t) this->offset_operation_[1];
	const size_t kk = (size_t) k + (size_t) this->offset_operation_[2];

#ifdef DEBUG
	for( int n=0; n<this->mpi_topology_->nproc; n++ ) {
		if( n == this->mpi_topology_->rank ) {
			printf("%d: PField::indexOperationToLocal\n", n);
			printf("  i, j, k = %d, %d, %d\n", i,j,k);
			printf("  (size_t)i, (size_t)j, (size_t)k = %ld, %ld, %ld\n", (size_t)i,(size_t)j,(size_t)k);
			printf("  ii, jj, kk = %ld, %ld, %ld\n", ii,jj,kk);
		}
		MPI_Barrier(this->mpi_topology_->comm);
	}
#endif
	return ((size_t) this->dims_local_[1] * ii + jj) * (size_t) this->dims_local_[2] + kk;
}

void PField::setSynchronized( bool synchronized )
{
	this->synchronized_ = synchronized;
}

// Set the grid pointers for field
/********************************************************************/
void PField::setGridLocal(const double *x_local, const double *y_local, const double *z_local)
/********************************************************************/
{
	this->x_local_ = x_local;
	this->y_local_ = y_local;
	this->z_local_ = z_local;
}

/********************************************************************/
void PField::setDataLocal( const size_t &index, const double &a )
/********************************************************************/
{
	// Set the single elemenetn
	this->data_local[index] = a;
	// Set unsynchronized
	this->synchronized_ = false;
}

/********************************************************************/
void PField::setDataOperation( const float *a)
/********************************************************************/
{
	*this = a;
	// Set unsynchronized
	this->synchronized_ = false;
}

/********************************************************************/
void PField::setDataOperation( const double *a)
/********************************************************************/
{
	*this = a;
	// Set unsynchronized
	this->synchronized_ = false;
}

/********************************************************************/
void PField::setDataOperation( const int &i, const int &j, const int &k, const double &a )
/********************************************************************/
{
	// Set the single elemenetn
	this->data_local[this->indexOperationToLocal(i, j, k)] = a;
	// Set unsynchronized
	this->synchronized_ = false;
}

//////////////////////////////////////////////////////////////////////
/// OPERATORS
//////////////////////////////////////////////////////////////////////

/********************************************************************/
PField& PField::operator=(const PField &a)
/********************************************************************/
{
	// Make sure they don't point to the same address
	if (this != &a) {
		// Copy the field data
		memcpy(this->data_local, a.data_local,
		       sizeof(*this->data_local) * this->getSizeLocal());
	}
	// Set unsynchronized
	this->synchronized_ = false;
	return *this;
}

/*
 * Sets the data on the operation domain. The size of "a" must be of
 * size getSizeOperation().
 */ 
/********************************************************************/
PField &PField::operator=( const double *a)
/********************************************************************/
{
	size_t index_operation;

	index_operation=0;
	PFIELD_LOOP_OPERATION_TO_LOCAL(this)
		this->data_local[_index] = a[index_operation++];
	PFIELD_LOOP_END

	// Set unsynchronized
	this->synchronized_ = false;
	return *this;
}

/*
 * Sets the data on the operation domain. The size of "a" must be of
 * size getSizeOperation().
 */ 
/********************************************************************/
PField &PField::operator=( const float *a)
/********************************************************************/
{
	size_t index_operation;

	index_operation=0;
	PFIELD_LOOP_OPERATION_TO_LOCAL(this)
		this->data_local[_index] = (double)a[index_operation++];
	PFIELD_LOOP_END

	// Set unsynchronized
	this->synchronized_ = false;
	return *this;
}


/********************************************************************/
PField& PField::operator=(float c)
/********************************************************************/
{
	*this = (double)c;
	return *this;
}

/********************************************************************/
PField& PField::operator=(double c)
/********************************************************************/
{
	std::fill_n( this->data_local, this->getSizeLocal(), c );

	// Set unsynchronized
	this->synchronized_ = false;
	return *this;
}

/********************************************************************/
PField& PField::operator+=(const PField &a)
/********************************************************************/
{

	PFIELD_LOOP_OPERATION_TO_LOCAL(this)
		this->data_local[_index] += a.data_local[_index];
	PFIELD_LOOP_END

	// Set unsynchronized
	this->synchronized_ = false;
	return *this;
}

/*
 * Sets the data on the operation domain. The size of "a" must be of
 * size getSizeOperation().
 */ 
/********************************************************************/
PField &PField::operator+=( const double *a)
/********************************************************************/
{
	size_t index_operation;

	index_operation=0;
	PFIELD_LOOP_OPERATION_TO_LOCAL(this)
		this->data_local[_index] += a[index_operation++];
	PFIELD_LOOP_END

	// Set unsynchronized
	this->synchronized_ = false;
	return *this;
}

/*
 * Sets the data on the operation domain. The size of "a" must be of
 * size getSizeOperation().
 */ 
/********************************************************************/
PField &PField::operator+=( const float *a)
/********************************************************************/
{
	size_t index_operation;

	index_operation=0;
	PFIELD_LOOP_OPERATION_TO_LOCAL(this)
		this->data_local[_index] += (double)a[index_operation++];
	PFIELD_LOOP_END

	// Set unsynchronized
	this->synchronized_ = false;
	return *this;
}
/********************************************************************/
PField& PField::operator+=(double c)
/********************************************************************/
{

	PFIELD_LOOP_OPERATION_TO_LOCAL(this)
		this->data_local[_index] += c;
	PFIELD_LOOP_END

	// Set unsynchronized
	this->synchronized_ = false;

	return *this;
}

/********************************************************************/
PField& PField::operator-=(const PField &a)
/********************************************************************/
{

	PFIELD_LOOP_OPERATION_TO_LOCAL(this)
		this->data_local[_index] -= a.data_local[_index];
	PFIELD_LOOP_END

	// Set unsynchronized
	this->synchronized_ = false;
	return *this;
}

/*
 * Sets the data on the operation domain. The size of "a"
 * must be of size getSizeOperation().
 */ 
/********************************************************************/
PField &PField::operator-=( const double *a)
/********************************************************************/
{
	size_t index_operation;

	index_operation=0;
	PFIELD_LOOP_OPERATION_TO_LOCAL(this)
		this->data_local[_index] -= a[index_operation++];
	PFIELD_LOOP_END

	// Set unsynchronized
	this->synchronized_ = false;
	return *this;
}

/*
 * Sets the data on the operation domain. The size of "a"
 * must be of size getSizeOperation().
 */ 
/********************************************************************/
PField &PField::operator-=( const float *a)
/********************************************************************/
{
	size_t index_operation;

	index_operation=0;
	PFIELD_LOOP_OPERATION_TO_LOCAL(this)
		this->data_local[_index] -= (double)a[index_operation++];
	PFIELD_LOOP_END

	// Set unsynchronized
	this->synchronized_ = false;
	return *this;
}

/********************************************************************/
PField& PField::operator-=(double c)
/********************************************************************/
{

	PFIELD_LOOP_OPERATION_TO_LOCAL(this)
		this->data_local[_index] -= c;
	PFIELD_LOOP_END

	// Set unsynchronized
	this->synchronized_ = false;

	return *this;
}

/********************************************************************/
PField& PField::operator*=(const PField &a)
/********************************************************************/
{
	PFIELD_LOOP_OPERATION_TO_LOCAL(this)
		this->data_local[_index] *= a.data_local[_index];
	PFIELD_LOOP_END

	// Set unsynchronized
	this->synchronized_ = false;

	return *this;
}

/*
 * Sets the data on the operation domain. The size of "a"
 * must be of size getSizeOperation().
 */ 
/********************************************************************/
PField &PField::operator*=( const double *a)
/********************************************************************/
{
	size_t index_operation;

	index_operation=0;
	PFIELD_LOOP_OPERATION_TO_LOCAL(this)
		this->data_local[_index] *= a[index_operation++];
	PFIELD_LOOP_END

	// Set unsynchronized
	this->synchronized_ = false;
	return *this;
}

/*
 * Sets the data on the operation domain. The size of "a"
 * must be of size getSizeOperation().
 */ 
/********************************************************************/
PField &PField::operator*=( const float *a)
/********************************************************************/
{
	size_t index_operation;

	index_operation=0;
	PFIELD_LOOP_OPERATION_TO_LOCAL(this)
		this->data_local[_index] *= (double)a[index_operation++];
	PFIELD_LOOP_END

	// Set unsynchronized
	this->synchronized_ = false;
	return *this;
}

/********************************************************************/
PField& PField::operator*=(float c)
/********************************************************************/
{
	*this *= (double)c;
	return *this;
}
/********************************************************************/
PField& PField::operator*=(double c)
/********************************************************************/
{

	PFIELD_LOOP_OPERATION_TO_LOCAL(this)
		this->data_local[_index] *= c;
	PFIELD_LOOP_END

	// Set unsynchronized
	this->synchronized_ = false;

	return *this;
}

/********************************************************************/
PField& PField::operator/=(const PField &a)
/********************************************************************/
{

	PFIELD_LOOP_OPERATION_TO_LOCAL(this)
		this->data_local[_index] /= a.data_local[_index];
	PFIELD_LOOP_END

	// Set unsynchronized
	this->synchronized_ = false;

	return *this;
}

/*
 * Sets the data on the operation domain. The size of "a"
 * must be of size getSizeOperation().
 */ 
/********************************************************************/
PField &PField::operator/=( const double *a)
/********************************************************************/
{
	size_t index_operation;

	index_operation=0;
	PFIELD_LOOP_OPERATION_TO_LOCAL(this)
		this->data_local[_index] /= a[index_operation++];
	PFIELD_LOOP_END

	// Set unsynchronized
	this->synchronized_ = false;
	return *this;
}

/*
 * Sets the data on the operation domain. The size of "a"
 * must be of size getSizeOperation().
 */ 
/********************************************************************/
PField &PField::operator/=( const float *a)
/********************************************************************/
{
	size_t index_operation;

	index_operation=0;
	PFIELD_LOOP_OPERATION_TO_LOCAL(this)
		this->data_local[_index] /= (double)a[index_operation++];
	PFIELD_LOOP_END

	// Set unsynchronized
	this->synchronized_ = false;
	return *this;
}

/********************************************************************/
PField& PField::operator/=(double c)
/********************************************************************/
{

	PFIELD_LOOP_OPERATION_TO_LOCAL(this)
		this->data_local[_index] /= c;
	PFIELD_LOOP_END

	// Set unsynchronized
	this->synchronized_ = false;

	return *this;
}

//////////////////////////////////////////////////////////////////////
/// MEMBER FUNCTIONS
//////////////////////////////////////////////////////////////////////

/********************************************************************/
PField &PField::add(const PField &a, const PField &b)
/********************************************************************/
{

#ifdef BOUNDS_CHECK
	size_t N = this->getSizeOperation();
	// First check that the fields are the same size
	if ( N != a.getSizeOperation() || N != b.getSizeOperation() ) {
		cout << "Mismatch in field sizes" << endl;
		exit (EXIT_FAILURE);
	}
#endif

	PFIELD_LOOP_OPERATION_TO_LOCAL(this)
		this->data_local[_index] = a.data_local[_index] + b.data_local[_index];
	PFIELD_LOOP_END

	// Set unsynchronized
	this->synchronized_ = false;
	return *this;
}

/********************************************************************/
PField &PField::sub(const PField &a, const PField &b)
/********************************************************************/
{

#ifdef BOUNDS_CHECK
	size_t N = this->getSizeOperation();
	// First check that the fields are the same size
	if ( N != a.getSizeOperation() || N != b.getSizeOperation() ) {
		cout << "Mismatch in field sizes" << endl;
		exit (EXIT_FAILURE);
	}
#endif

	PFIELD_LOOP_OPERATION_TO_LOCAL(this)
		this->data_local[_index] = a.data_local[_index] - b.data_local[_index];
	PFIELD_LOOP_END

	// Set unsynchronized
	this->synchronized_ = false;
	return *this;
}

/********************************************************************/
PField &PField::mul(const PField &a, const PField &b)
/********************************************************************/
{

#ifdef BOUNDS_CHECK
	size_t N = this->getSizeOperation();
	// First check that the fields are the same size
	if ( N != a.getSizeOperation() || N != b.getSizeOperation() ) {
		cout << "Mismatch in field sizes" << endl;
		exit (EXIT_FAILURE);
	}
#endif

	PFIELD_LOOP_OPERATION_TO_LOCAL(this)
		this->data_local[_index] = a.data_local[_index] * b.data_local[_index];
	PFIELD_LOOP_END

	// Set unsynchronized
	this->synchronized_ = false;
	return *this;
}

/********************************************************************/
PField &PField::div(const PField &a, const PField &b)
/********************************************************************/
{

#ifdef BOUNDS_CHECK
	size_t N = this->getSizeOperation();
	// First check that the fields are the same size
	if ( N != a.getSizeOperation() || N != b.getSizeOperation() ) {
		cout << "Mismatch in field sizes" << endl;
		exit (EXIT_FAILURE);
	}
#endif

	PFIELD_LOOP_OPERATION_TO_LOCAL(this)
		this->data_local[_index] = a.data_local[_index] / b.data_local[_index];
	PFIELD_LOOP_END

	// Set unsynchronized
	this->synchronized_ = false;
	return *this;
}

/********************************************************************/
PField &PField::sqrt(const PField &a)
/********************************************************************/
{

#ifdef BOUNDS_CHECK
	size_t N = this->getSizeOperation();
	// First check that the fields are the same size
	if ( N != a.getSizeOperation() ) {
		cout << "Mismatch in field sizes" << endl;
		exit (EXIT_FAILURE);
	}
#endif

	PFIELD_LOOP_OPERATION_TO_LOCAL(this)
		this->data_local[_index] = std::sqrt(a.data_local[_index]);
	PFIELD_LOOP_END

	// Set unsynchronized
	this->synchronized_ = false;
	return *this;
}

//////////////////////////////////////////////////////////////////////
/// DERIVATIVES (PUBLIC)
//////////////////////////////////////////////////////////////////////

// Public derivative functions use function pointers which call the private
// worker derivative functions. Note that mixed derivatives are computed by
// first computing the highest uni-directional (pure) derivative. Compound
// operations are then used to compute the subsequent derivatives. Hence we
// are not using pure mixed derivatives here.

/********************************************************************/
PField &PField::ddx(PField &a)
/********************************************************************/
{

#ifdef VERBOSE
	printf("%d: PField::ddx: entering\n", this->mpi_topology_->rank);
#endif

	// First check that the finite difference class has been initialized.
	int ierr=0;
	if( this->finite_diff_ == NULL ) {
		printf("%d: PField::ddx: finite difference class not initialized. Must call finiteDiffInit before calling this function\n", this->mpi_topology_->rank);
		MPI_Abort(this->mpi_topology_->comm,ierr);
	}

	void (FiniteDiff::*dd)(int, int, const double *, int, double *);
	dd = &FiniteDiff::ddx;
	this->dndxn(dd, a);

#ifdef VERBOSE
	printf("%d: PField::ddx: exiting\n", this->mpi_topology_->rank);
#endif
	return *this;
}

/********************************************************************/
PField &PField::ddy(PField &a)
/********************************************************************/
{

#ifdef VERBOSE
	printf("%d: PField::ddy: entering\n", this->mpi_topology_->rank);
#endif

	// First check that the finite difference class has been initialized.
	int ierr=0;
	if( this->finite_diff_ == NULL ) {
		printf("%d: PField::ddy: finite difference class not initialized. Must call finiteDiffInit before calling this function\n", this->mpi_topology_->rank);
		MPI_Abort(this->mpi_topology_->comm,ierr);
	}

	void (FiniteDiff::*dd)(int, int, const double *, int, double *);
	dd = &FiniteDiff::ddy;
	this->dndyn(dd, a);

#ifdef VERBOSE
	printf("%d: PField::ddy: exiting\n", this->mpi_topology_->rank);
#endif
	return *this;
}

/********************************************************************/
PField &PField::ddz(PField &a)
/********************************************************************/
{

	// First check that the finite difference class has been initialized.
	int ierr=0;
	if( this->finite_diff_ == NULL ) {
		printf("%d: PField::ddz: finite difference class not initialized. Must call finiteDiffInit before calling this function\n", this->mpi_topology_->rank);
		MPI_Abort(this->mpi_topology_->comm,ierr);
	}

	void (FiniteDiff::*dd)(int, int, const double *, int, double *);
	dd = &FiniteDiff::ddz;
	this->dndzn(dd, a);

	return *this;
}

/********************************************************************/
PField &PField::d2dx2(PField &a)
/********************************************************************/
{

	// First check that the finite difference class has been initialized.
	int ierr=0;
	if( this->finite_diff_ == NULL ) {
		printf("%d: PField::d2dx2: finite difference class not initialized. Must call finiteDiffInit before calling this function\n", this->mpi_topology_->rank);
		MPI_Abort(this->mpi_topology_->comm,ierr);
	}

	void (FiniteDiff::*dd)(int, int, const double *, int, double *);
	dd = &FiniteDiff::d2dx2;
	this->dndxn(dd, a);

	return *this;
}
/********************************************************************/
PField &PField::d2dy2(PField &a)
/********************************************************************/
{

	// First check that the finite difference class has been initialized.
	int ierr=0;
	if( this->finite_diff_ == NULL ) {
		printf("%d: PField::d2dy2: finite difference class not initialized. Must call finiteDiffInit before calling this function\n", this->mpi_topology_->rank);
		MPI_Abort(this->mpi_topology_->comm,ierr);
	}

	void (FiniteDiff::*dd)(int, int, const double *, int, double *);
	dd = &FiniteDiff::d2dy2;
	this->dndyn(dd, a);

	return *this;
}

/********************************************************************/
PField &PField::d2dz2(PField &a)
/********************************************************************/
{
	// First check that the finite difference class has been initialized.
	int ierr=0;
	if( this->finite_diff_ == NULL ) {
		printf("%d: PField::d2dz2: finite difference class not initialized. Must call finiteDiffInit before calling this function\n", this->mpi_topology_->rank);
		MPI_Abort(this->mpi_topology_->comm,ierr);
	}

	void (FiniteDiff::*dd)(int, int, const double *, int, double *);
	dd = &FiniteDiff::d2dz2;
	this->dndzn(dd, a);

	return *this;
}

/********************************************************************/
PField &PField::d2dxy(PField &a)
/********************************************************************/
{
	// First check that the finite difference class has been initialized.
	int ierr=0;
	if( this->finite_diff_ == NULL ) {
		printf("%d: PField::d2dxy: finite difference class not initialized. Must call finiteDiffInit before calling this function\n", this->mpi_topology_->rank);
		MPI_Abort(this->mpi_topology_->comm,ierr);
	}

	void (FiniteDiff::*dd)(int, int, const double *, int, double *);
	dd = &FiniteDiff::ddx;
	this->dndxn(dd, a);
	dd = &FiniteDiff::ddy;
	this->dndyn(dd, *this);

	return *this;
}
/********************************************************************/
PField &PField::d2dxz(PField &a)
/********************************************************************/
{
	// First check that the finite difference class has been initialized.
	int ierr=0;
	if( this->finite_diff_ == NULL ) {
		printf("%d: PField::d2dxz: finite difference class not initialized. Must call finiteDiffInit before calling this function\n", this->mpi_topology_->rank);
		MPI_Abort(this->mpi_topology_->comm,ierr);
	}

	void (FiniteDiff::*dd)(int, int, const double *, int, double *);
	dd = &FiniteDiff::ddx;
	this->dndxn(dd, a);
	dd = &FiniteDiff::ddz;
	this->dndzn(dd, *this);

	return *this;
}
/********************************************************************/
PField &PField::d2dyz(PField &a)
/********************************************************************/
{
	// First check that the finite difference class has been initialized.
	int ierr=0;
	if( this->finite_diff_ == NULL ) {
		printf("%d: PField::d2dyz: finite difference class not initialized. Must call finiteDiffInit before calling this function\n", this->mpi_topology_->rank);
		MPI_Abort(this->mpi_topology_->comm,ierr);
	}

	void (FiniteDiff::*dd)(int, int, const double *, int, double *);
	dd = &FiniteDiff::ddy;
	this->dndyn(dd, a);
	dd = &FiniteDiff::ddz;
	this->dndzn(dd, *this);

	return *this;
}

//////////////////////////////////////////////////////////////////////
/// DERIVATIVES (PRIVATE)
//////////////////////////////////////////////////////////////////////

/********************************************************************/
void PField::dndxn(void (FiniteDiff::*dd)(int, int, const double *, int, double *), PField &a)
/********************************************************************/
{

#ifdef VERBOSE
	printf("%d: PField::dndxn: entering\n", this->mpi_topology_->rank);
#endif

	// Create buffer to store x-data
	const int nx_local = this->dims_local_[0];
	const int nx_operation = this->dims_operation_[0];
	const int ny = this->dims_operation_[1];
	const int nz = this->dims_operation_[2];

	// Offsets
	const int x_offset = this->offset_operation_[0];
	const int y_offset = this->offset_operation_[1];
	const int z_offset = this->offset_operation_[2];

#ifdef BOUNDS_CHECK
	// First check that the fields are the same size
	// if ( nx != a.dims_operation_[0] || ny != a.dims_operation_[1] || nz != a.dims_operation_[2] ) {
	//   cout << "Mismatch in field sizes" << endl;
	//   exit (EXIT_FAILURE);
	// }
#endif

	double *ax = new double[nx_local];
	double *dax = new double[nx_operation];

	// Make sure the fields are synchronized; the state should be the same across all processes
	if (a.getSynchronized() == false)
		a.synchronize();

	for (int j = 0; j < ny; j++) {
		for (int k = 0; k < nz; k++) {
			// Pack buffer; need the entire buffer; have to manually apply operation offsets
			for (int i = 0; i < nx_local; i++)
				ax[i] = a.data_local[a.indexLocal(i, j + y_offset, k + z_offset)];
			// Compute the derivative using the FiniteDiff class
#ifdef DEBUG
			for(int n=0; n<this->mpi_topology_->nproc; n++) {
				if( n == this->mpi_topology_->rank ) {
					printf("%d:PField::dndxn\n", n);
					printf("    nx_local     = %d\n", nx_local);
					printf("    nx_operation = %d\n", nx_operation);
					printf("    ny           = %d\n", ny);
					printf("    nz           = %d\n", nz);
					for (int i = 0; i < nx_local; i++) {
						printf("    ax[%d] = %lf\n", i, ax[i]);
					}
				}
				MPI_Barrier(this->mpi_topology_->comm);
			}
#endif
			(this->finite_diff_->*dd)(x_offset, nx_local, ax, nx_operation,	dax);
			// Unpack buffer
			for (int i = 0; i < nx_operation; i++)
				this->data_local[this->indexOperationToLocal(i, j, k)] = dax[i];

		}
	}

	// Synchronize this field
	this->synchronize();

	delete[] ax;
	delete[] dax;

#ifdef VERBOSE
	printf("%d: PField::dndxn: exiting\n", this->mpi_topology_->rank);
#endif

}

/********************************************************************/
void PField::dndyn(void (FiniteDiff::*dd)(int, int, const double *, int, double *), PField &a)
/********************************************************************/
{

#ifdef VERBOSE
	printf("%d: PField::dndyn: entering\n", this->mpi_topology_->rank);
#endif

	// Create buffer to store x-data
	const int nx = this->dims_operation_[0];
	const int ny_local = this->dims_local_[1];
	const int ny_operation = this->dims_operation_[1];
	const int nz = this->dims_operation_[2];

	// Offsets
	const int x_offset = this->offset_operation_[0];
	const int y_offset = this->offset_operation_[1];
	const int z_offset = this->offset_operation_[2];

#ifdef BOUNDS_CHECK
	// First check that the fields are the same size
	// if ( nx != a.dims_operation_[0] || ny != a.dims_operation_[1] || nz != a.dims_operation_[2] ) {
	//   cout << "Mismatch in field sizes" << endl;
	//   exit (EXIT_FAILURE);
	// }
#endif

	double *ay = new double[ny_local];
	double *day = new double[ny_operation];

	// Make sure the fields are synchronized; the state should be the same across all processes
	if (a.getSynchronized() == false)
		a.synchronize();

	for (int i = 0; i < nx; i++) {
		for (int k = 0; k < nz; k++) {
			// Pack buffer
			for (int j = 0; j < ny_local; j++)
				ay[j] = a.data_local[a.indexLocal(i + x_offset, j, k + z_offset)];
			// Compute the derivative using the FiniteDiff class
			(this->finite_diff_->*dd)(y_offset, ny_local, ay, ny_operation,	day);
			// Unpack buffer
			for (int j = 0; j < ny_operation; j++)
				this->data_local[this->indexOperationToLocal(i, j, k)] = day[j];
		}
	}

	// Synchronize this field
	this->synchronize();

	delete[] ay;
	delete[] day;

#ifdef VERBOSE
	printf("%d: PField::dndyn: exiting\n", this->mpi_topology_->rank);
#endif

}

/********************************************************************/
void PField::dndzn(void (FiniteDiff::*dd)(int, int, const double *, int, double *), PField &a)
/********************************************************************/
{

#ifdef VERBOSE
	printf("%d: PField::dndzn: entering\n", this->mpi_topology_->rank);
#endif

	// Create buffer to store x-data
	const int nx = this->dims_operation_[0];
	const int ny = this->dims_operation_[1];
	const int nz_local = this->dims_local_[2];
	const int nz_operation = this->dims_operation_[2];

	// Offsets
	const int x_offset = this->offset_operation_[0];
	const int y_offset = this->offset_operation_[1];
	const int z_offset = this->offset_operation_[2];

#ifdef BOUNDS_CHECK
	// First check that the fields are the same size
	// if ( nx != a.dims_operation_[0] || ny != a.dims_operation_[1] || nz != a.dims_operation_[2] ) {
	//   cout << "Mismatch in field sizes" << endl;
	//   exit (EXIT_FAILURE);
	// }
#endif

	double *az = new double[nz_local];
	double *daz = new double[nz_operation];

	// Make sure the fields are synchronized; the state should be the same across all processes
	if (a.getSynchronized() == false)
		a.synchronize();

	for (int i = 0; i < nx; i++) {
		for (int j = 0; j < ny; j++) {
			// Pack buffer
			for (int k = 0; k < nz_local; k++)
				az[k] = a.data_local[a.indexLocal(i + x_offset, j + y_offset, k)];
			// Compute the derivative using the FiniteDiff class
			(this->finite_diff_->*dd)(z_offset, nz_local, az, nz_operation,	daz);
			// Unpack buffer
			for (int k = 0; k < nz_operation; k++)
				this->data_local[this->indexOperationToLocal(i, j, k)] = daz[k];
		}
	}

	// Synchronize this field
	this->synchronize();

	delete[] az;
	delete[] daz;

#ifdef VERBOSE
	printf("%d: PField::dndzn: exiting\n", this->mpi_topology_->rank);
#endif

}

// 
// Procedures used for filtering the field
// 
//////////////////////////////////////////////////////////////////////
/// FILTERING (PUBLIC)
//////////////////////////////////////////////////////////////////////

PField &PField::filter( const int &filter_width )
{
	
	int ierr=0;

	// First make sure the filter width is an even value
	if( filter_width % 2 != 0 ) {
		printf("%d: PField::filter: filter width must be an even, positive integer\n", this->mpi_topology_->rank);
		MPI_Abort(this->mpi_topology_->comm,ierr);
	}	    
	// Also make sure that the filter width is not too wide since overlap data is needed
	if( filter_width/2 > this->rind_size_ ) {
		printf("%d: PField::filter: rind size too small to support filter width\n", this->mpi_topology_->rank);
		MPI_Abort(this->mpi_topology_->comm,ierr);
	}

	const int filter_width_half = filter_width / 2;

	// Determine which function to used to select the bounds for the filtering
	const int (PField::*x_filter_width_left)( const int &, const int &, const int &);
	const int (PField::*x_filter_width_right)( const int &, const int &, const int &);

	if( this->hasRind(0,-1) ) {
		x_filter_width_left = &PField::filter_width_left_rind;
	} else {
		x_filter_width_left = &PField::filter_width_left_boundary;
	}
	if( this->hasRind(0,1) ) {
		x_filter_width_right = &PField::filter_width_right_rind;
	} else {
		x_filter_width_right = &PField::filter_width_right_boundary;
	}

	const int (PField::*y_filter_width_left)( const int &, const int &, const int &);
	const int (PField::*y_filter_width_right)( const int &, const int &, const int &);

	if( this->hasRind(1,-1) ) {
		y_filter_width_left = &PField::filter_width_left_rind;
	} else {
		y_filter_width_left = &PField::filter_width_left_boundary;
	}
	if( this->hasRind(1,1) ) {
		y_filter_width_right = &PField::filter_width_right_rind;
	} else {
		y_filter_width_right = &PField::filter_width_right_boundary;
	}

	const int (PField::*z_filter_width_left)( const int &, const int &, const int &);
	const int (PField::*z_filter_width_right)( const int &, const int &, const int &);

	if( this->hasRind(2,-1) ) {
		z_filter_width_left = &PField::filter_width_left_rind;
	} else {
		z_filter_width_left = &PField::filter_width_left_boundary;
	}
	if( this->hasRind(2,1) ) {
		z_filter_width_right = &PField::filter_width_right_rind;
	} else {
		z_filter_width_right = &PField::filter_width_right_boundary;
	}

	static const double eight_inv = 1.0L / 8.0L;
	
	// Need a buffer to store the filtered operation domain data
	double *data_domain = new double[this->getSizeOperation()];

	// Now synchronize the field
	this->synchronize();

	// Need to loop over the entire operation domain
	PFIELD_LOOP_OPERATION(this)

		// Obtain the widths of the filter stencil for the left and right parts
		const int qx_left  = (this->*x_filter_width_left) (0, _i, filter_width_half);
	        const int qx_right = (this->*x_filter_width_right)(0, _i, filter_width_half);
		const int qy_left  = (this->*y_filter_width_left) (1, _j, filter_width_half);
		const int qy_right = (this->*y_filter_width_right)(1, _j, filter_width_half);
		const int qz_left  = (this->*z_filter_width_left) (2, _k, filter_width_half);
		const int qz_right = (this->*z_filter_width_right)(2, _k, filter_width_half);

		// Peform summation over the filter stencil
		data_domain[_index] = 0.0L;
		for(int i= _i-qx_left; i< _i + qx_right; i++ ) {

			const int i_local = _i + this->offset_operation_[0];
			const double dx = this->x_local_[i_local+1] - this->x_local_[i_local];

			for(int j= _j-qy_left; j< _j + qy_right; j++ ) {

				const int j_local = _j + this->offset_operation_[1];
				const double dy = this->y_local_[j_local+1] - this->y_local_[j_local];

				for(int k= _k-qz_left; k< _k + qz_right; k++ ) {

					const int k_local = _k + this->offset_operation_[2];
					const double dz = this->z_local_[k_local+1] - this->z_local_[k_local];

					data_domain[_index] += 
						( this->data_local[this->indexLocal(i_local  ,j_local  ,k_local  )] + 
						  this->data_local[this->indexLocal(i_local  ,j_local  ,k_local+1)] + 
						  this->data_local[this->indexLocal(i_local  ,j_local+1,k_local  )] + 
						  this->data_local[this->indexLocal(i_local  ,j_local+1,k_local+1)] + 
						  this->data_local[this->indexLocal(i_local+1,j_local  ,k_local  )] + 
						  this->data_local[this->indexLocal(i_local+1,j_local  ,k_local+1)] + 
						  this->data_local[this->indexLocal(i_local+1,j_local+1,k_local  )] + 
						  this->data_local[this->indexLocal(i_local+1,j_local+1,k_local+1)] 
						  ) * dx * dy * dz;
				}
			}
		}

		const double x_filter_left  = this->x_local_[_i - qx_left  + this->offset_operation_[0] ];
		const double x_filter_right = this->x_local_[_i + qx_right + this->offset_operation_[0] ];
		const double y_filter_left  = this->y_local_[_j - qy_left  + this->offset_operation_[1] ];
		const double y_filter_right = this->y_local_[_j + qy_right + this->offset_operation_[1] ];
		const double z_filter_left  = this->z_local_[_k - qz_left  + this->offset_operation_[2] ];
		const double z_filter_right = this->z_local_[_k + qz_right + this->offset_operation_[2] ];

		// Now apply the normalization
		data_domain[_index] *= eight_inv / 
			               ( x_filter_right - x_filter_left ) / 
			               ( y_filter_right - y_filter_left ) / 
			               ( z_filter_right - z_filter_left ); 


	PFIELD_LOOP_END

	// Need to loop over the entire operation domain to set data_local
	size_t index_operation=0;
	PFIELD_LOOP_OPERATION_TO_LOCAL(this)
		this->data_local[_index] = data_domain[index_operation++];
	PFIELD_LOOP_END

	delete [] data_domain;

	return *this;
}


//
// The computation procedures compute and return information regarding data
// members of the class
//
////////////////////////////////////////////////////////////////////////////////
/// COMPUTATIONS
////////////////////////////////////////////////////////////////////////////////

//
// Computes the dimensions/size of the MPI topology for a given number of
// processes and domain decomposition dimensions. The size of
// mpi_topology_dims matches PFIELD_NDIMS.
//
// Inputs:
//   nproc - number of total MPI processes
//   mpi_decomp_ndims - number of decomposition dimensions (1 for slab decomposition, 2 for pencil, etc.)
// Returns:
//   mpi_topology_dims[PFIELD_NDIMS] - topology dimensions in each field direction
//
/********************************************************************/
int *PField::computeMPITopologyDims(int nproc, int mpi_decomp_ndims) const
/********************************************************************/
{

	// initialize the MPI topology dimensions to 1
	int *mpi_topology_dims = new int[PFIELD_NDIMS];
	std::fill_n(mpi_topology_dims, PFIELD_NDIMS, 1);

	if (mpi_decomp_ndims == 1) {

		mpi_topology_dims[0] = nproc;

	} else if (mpi_decomp_ndims == 2) {

		int M = std::min( (int) std::ceil( std::sqrt(nproc * (double) this->dims_[0] / this->dims_[1])), nproc);
		int N = nproc / M;

		while (1) {
			if (N * M == nproc) {
				break;
			} else if (N * M < nproc) {
				M++;
			} else {
				std::cout << "Unable to set topology" << std::endl;
				exit(EXIT_FAILURE);
			}
			N = nproc / M;
		}

		mpi_topology_dims[0] = M;
		mpi_topology_dims[1] = N;

	} else {

		std::cout << "Currently only MPI topologies up to 2 dimensions are supported" << std::endl;
		exit(EXIT_FAILURE);

	}

	return mpi_topology_dims;

}

//
// The assigner procedures compute and assign data members of the class
//
////////////////////////////////////////////////////////////////////////////////
/// ASSIGNERS
////////////////////////////////////////////////////////////////////////////////

/********************************************************************/
void PField::assignMPITopology()
/********************************************************************/
{

	int mpi_decomp_ndims;

	// Determine the MPI coordinate dimensions
	if (this->field_decomp_ == FIELD_DECOMP_SLAB) {
		mpi_decomp_ndims = 1;
	} else if (this->field_decomp_ == FIELD_DECOMP_PENCIL) {
		mpi_decomp_ndims = 2;
	} else {
		cout << "only slab or pencil decompositions currently supported" << endl;
		exit(EXIT_FAILURE);
	}

	int nproc;

	// Get the total number of procs
	MPI_Comm_size(MPI_COMM_WORLD, &nproc);

	int *mpi_topology_dims = this->computeMPITopologyDims(nproc, mpi_decomp_ndims);

	// Create the MPI topology struct; this contains all the comm, rank, nproc, neighbors, etc.
	this->mpi_topology_ = MPITopologyNew(MPI_COMM_WORLD, PFIELD_NDIMS, mpi_topology_dims, this->periodic_);

	delete[] mpi_topology_dims;

#ifdef VERBOSE
	for(int n=0; n<this->mpi_topology_->nproc; n++) {
		if( n == this->mpi_topology_->rank ) {
			printf("%d: PField::assignMPITopology\n",this->mpi_topology_->rank);
			printf("    coords: %d, %d\n",this->mpi_topology_->coords[0], this->mpi_topology_->coords[1]);
			printf("    neighbor_prev: %d, %d\n",this->mpi_topology_->neighbor_prev[0], this->mpi_topology_->neighbor_prev[1]);
			printf("    neighbor_next: %d, %d\n",this->mpi_topology_->neighbor_next[0], this->mpi_topology_->neighbor_next[1]);
		}
		MPI_Barrier(this->mpi_topology_->comm);
	}
#endif

}

/********************************************************************/
void PField::assignDimsAndOffsets()
/********************************************************************/
{

	int *zero3 = new int[PFIELD_NDIMS]; // Zero vector of dimension PFIELD_NDIMS
	std::fill_n(zero3, PFIELD_NDIMS, 0);

	// Initialize the local/operation dimensions and offset
	memcpy(this->dims_local_, this->dims_, sizeof(*this->dims_) * PFIELD_NDIMS);
	memcpy(this->offset_local_, zero3, sizeof(*this->offset_local_) * PFIELD_NDIMS);

	memcpy(this->dims_operation_, this->dims_, sizeof(*this->dims_) * PFIELD_NDIMS);
	memcpy(this->offset_operation_, zero3, sizeof(*this->offset_local_) * PFIELD_NDIMS);

	delete[] zero3;

	// For each dimension in the MPI topology; determin what the dimensions and
	// offsets are of the local and operation domains
	for (int n = 0; n < this->mpi_topology_->ndims; n++) {

		int p = this->dims_[n] % this->mpi_topology_->dims[n]; // Remainder of points (surplus)

		int d = this->dims_[n] / this->mpi_topology_->dims[n];
		int s = (this->mpi_topology_->coords[n] < p) ? 1 : 0; // Take one of the surplus values if in low enough rank

		// Check that the rind size is not larger than the size of the smallest dimension
		if (d < this->rind_size_) {
			cout << "Domain decomposition too fine to support derivatives of order " << this->operator_order_ << endl;
			exit(EXIT_FAILURE);
		}

		this->dims_local_[n] = d + s; // Local domain size without rind data
		// Offset of the local domain with respect to the global domain
		this->offset_local_[n] = this->mpi_topology_->coords[n] * d + std::min(this->mpi_topology_->coords[n], p);

		// The operation domain is the local domain without rind data
		this->dims_operation_[n] = this->dims_local_[n];

		// Now determine the full local dimensions with the rind
		if (this->mpi_topology_->neighbor_next[n] >= 0)
			this->dims_local_[n] += this->rind_size_; // If we have a neighbor then add the rind points
		if (this->mpi_topology_->neighbor_prev[n] >= 0) {
			this->dims_local_[n] += this->rind_size_;
			this->offset_local_[n] -= this->rind_size_; // Move the local offset position back to compensate for the rind points
			this->offset_operation_[n] = this->rind_size_; // Set the operation offset distance relative to the local domain
		}
	}
}

////////////////////////////////////////////////////////////////////////////////
/// MPI TRANSFER FUNCTIONS
////////////////////////////////////////////////////////////////////////////////

/*
 * Procedure to check if the field has a rind for a given dimension and location
 */
/********************************************************************/
bool PField::hasRind(int dim, int location) const
/********************************************************************/
{
	if (location == -1) {
		return (this->mpi_topology_->neighbor_prev[dim] >= 0) ? true : false;
	} else if( location == 1 ) {
		return (this->mpi_topology_->neighbor_next[dim] >= 0) ? true : false;
	} else {
		printf("%d: PField::hasRind: incorrect location provided\n", this->mpi_topology_->rank);
		int ierr=0;
		MPI_Abort(this->mpi_topology_->comm, ierr);
		return false; // This shouldn't be reached but put here to avoid compiler warning
	}
}

/********************************************************************/
void PField::synchronize()
/********************************************************************/
{

#ifdef VERBOSE
        printf("%d: PField::synchronize\n",this->mpi_topology_->rank);
	printf("%d:     field_decomp = %d\n", this->mpi_topology_->rank, this->field_decomp_);
#endif
	if( this->synchronized_ == true ) return;

	// x pass
	if (this->field_decomp_ == FIELD_DECOMP_SLAB || 
	    this->field_decomp_ == FIELD_DECOMP_PENCIL) {
		this->synchronizeDimension(0);
	}

	if (this->field_decomp_ == FIELD_DECOMP_PENCIL) {
		this->synchronizeDimension(1);
	}

	// Set that the field is synchronized
	this->synchronized_ = true;

}

/********************************************************************/
void PField::synchronizeDimension(int dim)
/********************************************************************/
{
	MPI_Request prev_request=0, next_request=0;
	MPI_Status status;

	MPITopology_t *mpi_topology = this->mpi_topology_;

	const size_t prev_buffer_size = this->getSizeRind(dim, -1);
	const size_t next_buffer_size = this->getSizeRind(dim, 1);
	double *prev_buffer = this->createRindBuffer(dim, -1);
	double *next_buffer = this->createRindBuffer(dim, 1);

#ifdef VERBOSE
	printf("%d: PField::synchronizeDimension: has rind next = %d\n", mpi_topology->rank, this->hasRind(dim,1));
	printf("%d: PField::synchronizeDimension: has rind prev = %d\n", mpi_topology->rank, this->hasRind(dim,-1));                                           
#endif 
	/*
	 * First handshakes
	 */
	// Recv next
	if( this->hasRind( dim, 1 ) ) {
#ifdef VERBOSE
	printf("%d: PField::synchronizeDimension: receiving from next = %d\n", mpi_topology->rank, mpi_topology->neighbor_next[dim]);
#endif
		MPI_Irecv(next_buffer, next_buffer_size, MPI_DOUBLE,
			  mpi_topology->neighbor_next[dim], 1, mpi_topology->comm,
			  &next_request);
	}

	if( this->hasRind( dim, -1 ) ) {
#ifdef VERBOSE
		printf("%d: PField::synchronizeDimension: sending to prev = %d\n", mpi_topology->rank, mpi_topology->neighbor_prev[dim]);
#endif
		// Send prev
		this->packRindBuffer(dim, -1, prev_buffer);

		MPI_Isend(prev_buffer, prev_buffer_size, MPI_DOUBLE,
			  mpi_topology->neighbor_prev[dim], 1, mpi_topology->comm,
			  &prev_request);
	}

	if( this->hasRind( dim, 1 ) ) {
#ifdef VERBOSE
		printf("%d: PField::synchronizeDimension: waiting to receive from next\n", mpi_topology->rank);
#endif
		MPI_Wait(&next_request, &status);
		this->unpackRindBuffer(dim, 1, next_buffer);
	}

	if( this->hasRind( dim, -1 ) ) {
		MPI_Wait(&prev_request, &status);
	}

	/*
	 * Second handshakes
	 */
	if( this->hasRind(dim, -1) ) {
		// Recv prev
#ifdef VERBOSE
		printf("%d: PField::synchronizeDimension: receiving from prev = %d\n", mpi_topology->rank, mpi_topology->neighbor_prev[dim]);
#endif
		MPI_Irecv(prev_buffer, prev_buffer_size, MPI_DOUBLE,
			  mpi_topology->neighbor_prev[dim], 2, mpi_topology->comm,
			  &prev_request);
	}

	if( this->hasRind(dim, 1 ) ) {
#ifdef VERBOSE
		printf("%d: PField::synchronizeDimension: sending to next\n", mpi_topology->rank);
#endif
	// Send next
		this->packRindBuffer(dim, 1, next_buffer);
		MPI_Isend(next_buffer, next_buffer_size, MPI_DOUBLE,
			  mpi_topology->neighbor_next[dim], 2, mpi_topology->comm,
			  &next_request);
	}
	
	if( this->hasRind(dim, -1) ) {
#ifdef VERBOSE
		printf("%d: PField::synchronizeDimension: receiving from prev\n", mpi_topology->rank);
#endif
		MPI_Wait(&prev_request, &status);
		this->unpackRindBuffer(dim, -1, prev_buffer);
	}

	if( this->hasRind(dim, 1) ) {
		MPI_Wait(&next_request, &status);
	}

	if( prev_buffer != NULL ) delete[] prev_buffer;
	if( next_buffer != NULL ) delete[] next_buffer;
	
}

/********************************************************************/
double *PField::createRindBuffer(int dim, int location) const
/********************************************************************/
{
	if( this->hasRind(dim,location) ) {
		return new double[this->getSizeRind(dim, location)];
	} else {
		return NULL;
	}
}

//
// This procedure packs the data buffer for MPI transfers.
//
// Inputs:
//   dim - dimension to pack buffer for (0-x, 1-y, 2-z)
//   location - location in domain (-1-bottom/left, 1-top/right)
// Outputs:
//   rind_buffer - buffer to pack data into
//
// Note: rind_buffer must be the correct size. Use createRindBuffer to
// allocate the buffer.
//
/********************************************************************/
void PField::packRindBuffer(int dim, int location, double *rind_buffer) const
/********************************************************************/
{

#ifdef VERBOSE
	printf("%d: PField::packRindBuffer: entering\n", this->mpi_topology_->rank);
#endif
	int imin, imax;
	int jmin, jmax;
	int kmin, kmax;

	// If there is no rind do nothing
	if ( !this->hasRind(dim, location) ) {
		return;
	}

	imin = 0;
	imax = this->dims_operation_[0] - 1;

	jmin = 0;
	jmax = this->dims_operation_[1] - 1;
	
	kmin = 0;
	kmax = this->dims_operation_[2] - 1;

	if (dim == 0) {

		if (location == -1) {
			imin = 0;
		} else {
			imin = this->dims_operation_[0] - this->rind_size_;
		}
		imax = imin + this->rind_size_ - 1;

	} else if (dim == 1) {

		if (location == -1) {
			jmin = 0;
		} else {
			jmin = this->dims_operation_[1] - this->rind_size_;
		}
		jmax = jmin + this->rind_size_ - 1;

	} else if (dim == 2) {

		if (location == -1) {
			kmin = 0;
		} else {
			kmin = this->dims_operation_[2] - this->rind_size_;
		}
		kmax = kmin + this->rind_size_ - 1;

	} else {
		printf("%d: PField::packRindBuffer: unable to pack buffer: incorrect dimension specified\n", this->mpi_topology_->rank);
		int ierr=0;
		MPI_Abort(this->mpi_topology_->comm, ierr);
	}

	size_t index = 0;
	int i=0, j=0, k=0;
	for (i = imin; i <= imax; i++) {
		for (j = jmin; j <= jmax; j++) {
			for (k = kmin; k <= kmax; k++) {
				rind_buffer[index] =
					this->data_local[this->indexOperationToLocal(i, j, k)];
				index++;
			}
		}
	}

	// Perform sanity check
	if (index != this->getSizeRind(dim, location)) {
		printf("%d: PField::packRindBuffer: mismatch in expected rind buffer size\n", this->mpi_topology_->rank);
		int ierr=0;
		MPI_Abort(this->mpi_topology_->comm, ierr);
	}

#ifdef VERBOSE
	printf("%d: PField::packRindBuffer: exiting\n", this->mpi_topology_->rank);
#endif

}

//
// This procedure unpacks the data buffer for MPI transfers.
//
// Inputs:
//   dim - dimension to pack buffer for (0-x, 1-y, 2-z)
//   location - location in domain (-1-bottom/left, 1-top/right)
//   rind_buffer - buffer to pack data into
//
// Note: rind_buffer must be the correct size. Use createRindBuffer to
// allocate the buffer.
//
/********************************************************************/
void PField::unpackRindBuffer(int dim, int location, const double *rind_buffer)
/********************************************************************/
{

#ifdef VERBOSE
	printf("%d: PField::unpackRindBuffer: entering\n", this->mpi_topology_->rank);
#endif

	int istart, isize;
	int jstart, jsize;
	int kstart, ksize;

	// If there is no rind do nothing
	if ( !this->hasRind(dim, location) ) {
		return;
	}

	// Preset the starting indicies and sizes
	istart = this->offset_operation_[0];
	isize = this->dims_operation_[0];

	jstart = this->offset_operation_[1];
	jsize = this->dims_operation_[1];

	kstart = this->offset_operation_[2];
	ksize = this->dims_operation_[2];

	// Reset values asize_t dim direction
	if (dim == 0) {

		if (location == -1) {
			istart = 0;
		} else {
			istart = this->dims_local_[0] - this->rind_size_;
		}
		isize = this->rind_size_;

	} else if (dim == 1) {

		if (location == -1) {
			jstart = 0;
		} else {
			jstart = this->dims_local_[1] - this->rind_size_;
		}
		jsize = this->rind_size_;

	} else if (dim == 2) {

		if (location == -1) {
			kstart = 0;
		} else {
			kstart = this->dims_local_[2] - this->rind_size_;
		}
		ksize = this->rind_size_;

	} else {
		std::cout
			<< "PField::unpackRindBuffer: unable to pack buffer: incorrect dimension specified"
			<< std::endl;
		exit(EXIT_FAILURE);
	}

	size_t index = 0;
	for (int i = istart; i < istart + isize; i++) {
		for (int j = jstart; j < jstart + jsize; j++) {
			for (int k = kstart; k < kstart + ksize; k++) {
				this->data_local[this->indexLocal(i, j, k)] =
					rind_buffer[index++];
			}
		}
	}

	// Perform sanity check
	if (index != this->getSizeRind(dim, location)) {
		printf("%d: PField::packRindBuffer: mismatch in expected rind buffer size\n", this->mpi_topology_->rank);
		int ierr=0;
		MPI_Abort(this->mpi_topology_->comm, ierr);
	}

#ifdef VERBOSE
	printf("%d: PField::unpackRindBuffer: exiting\n", this->mpi_topology_->rank);
#endif

}

//////////////////////////////////////////////////////////////////////
/// NON-CLASS MEMBER FUNCTIONS
//////////////////////////////////////////////////////////////////////

}

