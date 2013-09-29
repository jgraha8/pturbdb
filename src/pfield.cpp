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

PField::PField(PField &g, bool copy_data_local)
{

	this->PFieldCopy( g, copy_data_local );

#ifdef VERBOSE
	printf("%d: PField::PField: constructed field from copy\n", this->getMPITopology()->rank);
#endif

}

PField::PField(PField &g)
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

       	delete[] this->data_local;
	delete this->finite_diff_;
	delete this->mpi_topology_;

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

void PField::PFieldCopy( PField &g, bool copy_data_local )
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

	// Copy pointers
	this->mpi_topology_   = g.getMPITopology();
	this->finite_diff_    = g.getFiniteDiff();

	this->x_local_        = g.getXLocal();
	this->y_local_        = g.getYLocal();
	this->z_local_        = g.getZLocal();

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
		printf("%d: PField::finiteDiffInit: class instance fd already created\n", this->mpi_topology_->rank);
		MPI_Abort(this->mpi_topology_->comm, ierr);
	}

	if (this->x_local_ == NULL || this->y_local_ == NULL
	    || this->z_local_ == NULL) {
		printf("%d: PField::finiteDiffInit: requires setGridLocal be called first\n", this->mpi_topology_->rank);
		MPI_Abort(this->mpi_topology_->comm, ierr);
	}

	// First initialize the fd class instance
	this->finite_diff_ = new FiniteDiff(this->dims_local_[0], this->x_local_,
					    this->dims_local_[1], this->y_local_, this->dims_local_[2],
					    this->z_local_, this->operator_order_);

#ifdef VERBOSE
	printf("%d: PField::finiteDiffInit: exiting\n", this->mpi_topology_->rank);
#endif

}

//////////////////////////////////////////////////////////////////////
/// GETTERS/SETTERS
//////////////////////////////////////////////////////////////////////

/********************************************************************/
long PField::getSize()
/********************************************************************/
{
	long N = 1;
	for (int n = 0; n < PFIELD_NDIMS; n++) {
		N *= this->dims_[n];
	}
	return N;
}

/********************************************************************/
long PField::getSizeLocal()
/********************************************************************/
{
	long N = 1;
	for (int n = 0; n < PFIELD_NDIMS; n++) {
		N *= this->dims_local_[n];
	}
	return N;
}

/********************************************************************/
long PField::getSizeOperation()
/********************************************************************/
{
	long N = 1;
	for (int n = 0; n < PFIELD_NDIMS; n++) {
		N *= this->dims_operation_[n];
	}
	return N;
}

/********************************************************************/
long PField::getSizeRind(int dim, int location)
/********************************************************************/
{
	// If the dimension and location does not have a rind, return 0
	if ( !this->hasRind(dim, location) ) 
		return 0;

	if (dim == 0) {
		return (long) this->rind_size_ * (long) this->dims_operation_[1]
			* (long) this->dims_operation_[2];
	} else if (dim == 1) {
		return (long) this->dims_operation_[0] * (long) this->rind_size_
			* (long) this->dims_operation_[2];
	} else if (dim == 2) {
		return (long) this->dims_operation_[0] * (long) this->dims_operation_[1]
			* (long) this->rind_size_;
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
FieldDecomp_t PField::getFieldDecomp()  { return this->field_decomp_;     };
MPITopology_t *PField::getMPITopology() { return this->mpi_topology_;     };
FiniteDiff *PField::getFiniteDiff()     { return this->finite_diff_;      };
double *PField::getXLocal()             { return this->x_local_;          };
double *PField::getYLocal()             { return this->y_local_;          };
double *PField::getZLocal()             { return this->z_local_;          };

int *PField::getFieldPeriodic()         { return this->periodic_;         };

int PField::getOperatorOrder()          { return this->operator_order_;   };

int *PField::getDims()                  { return this->dims_;             };
int *PField::getDimsLocal()             { return this->dims_local_;       };    
int *PField::getDimsOperation()         { return this->dims_operation_;   };

int *PField::getOffsetLocal()           { return this->offset_local_;     };
int *PField::getOffsetOperation()       { return this->offset_operation_; };

int PField::getRindSize()               { return this->rind_size_;        };

bool PField::getSynchronized()          { return this->synchronized_;     };

/*
 * Gets the x data on the operation grid by copying to a provided data
 * buffer. The size of "x" must be of size getDimsOperation().
 */ 
void PField::getXOperation( double *x)
{
	for (int i = 0; i < this->dims_operation_[0]; i++)
		x[i] = this->x_local_[i + this->offset_operation_[0]];
}
/*
 * Gets the y data on the operation grid by copying to a provided data
 * buffer. The size of "y" must be of size getDimsOperation().
 */ 
void PField::getYOperation( double *y)
{
	for (int j = 0; j < this->dims_operation_[1]; j++)
		y[j] = this->y_local_[j + this->offset_operation_[1]];
}
/*
 * Gets the z data on the operation grid by copying to a provided data
 * buffer. The size of "z" must be of size getDimsOperation().
 */ 
void PField::getZOperation( double *z)
{
	for (int k = 0; k < this->dims_operation_[2]; k++)
		z[k] = this->z_local_[k + this->offset_operation_[2]];
}


/*
 * Gets the data on the operation grid by copying to a provided data
 * buffer. The size of "a" must be of size getSizeOperation().
 */ 
void PField::getDataOperation( float *a)
{
	long index_local, index_operation;

	index_operation=0;
	for (int i = 0; i < this->dims_operation_[0]; i++) {
		for (int j = 0; j < this->dims_operation_[1]; j++) {
			for (int k = 0; k < this->dims_operation_[2]; k++) {
				index_local = this->indexOperationToLocal(i, j, k);
				a[index_operation] = (float)this->data_local[index_local];
				++index_operation;
			}
		}
	}
}

/*
 * Gets the data on the operation grid by copying to a provided data
 * buffer. The size of "a" must be of size getSizeOperation().
 */ 
void PField::getDataOperation( double *a)
{
	long index_local, index_operation;

	index_operation=0;
	for (int i = 0; i < this->dims_operation_[0]; i++) {
		for (int j = 0; j < this->dims_operation_[1]; j++) {
			for (int k = 0; k < this->dims_operation_[2]; k++) {
				index_local = this->indexOperationToLocal(i, j, k);
				a[index_operation] = this->data_local[index_local];
				++index_operation;
			}
		}
	}
}


// Array index for 3D indexing
/********************************************************************/
long PField::index(int i, int j, int k)
/********************************************************************/
{
	return ((long) this->dims_[1] * (long) i + (long) j) * (long) this->dims_[2]
		+ (long) k;
}

// Array index for 3D indexing
/********************************************************************/
long PField::indexLocal(int i, int j, int k)
/********************************************************************/
{
	return ((long) this->dims_local_[1] * (long) i + (long) j)
		* (long) this->dims_local_[2] + (long) k;
}

// Array index for 3D indexing
/********************************************************************/
long PField::indexOperation(int i, int j, int k)
/********************************************************************/
{
	return ((long) this->dims_operation_[1] * (long) i + (long) j)
		* (long) this->dims_operation_[2] + (long) k;
}

// Array index for 3D indexing
/********************************************************************/
long PField::indexOperationToLocal(int i, int j, int k)
/********************************************************************/
{
	const long ii = (long) i + (long) this->offset_operation_[0];
	const long jj = (long) j + (long) this->offset_operation_[1];
	const long kk = (long) k + (long) this->offset_operation_[2];

#ifdef DEBUG
	for( int n=0; n<this->mpi_topology_->nproc; n++ ) {
		if( n == this->mpi_topology_->rank ) {
			printf("%d: PField::indexOperationToLocal\n", n);
			printf("  i, j, k = %d, %d, %d\n", i,j,k);
			printf("  (long)i, (long)j, (long)k = %ld, %ld, %ld\n", (long)i,(long)j,(long)k);
			printf("  ii, jj, kk = %ld, %ld, %ld\n", ii,jj,kk);
		}
		MPI_Barrier(this->mpi_topology_->comm);
	}
#endif
	return ((long) this->dims_local_[1] * ii + jj) * (long) this->dims_local_[2]
		+ kk;
}

void PField::setSynchronized( bool synchronized ) 
{
	this->synchronized_ = synchronized;
}

// Set the grid pointers for field
/********************************************************************/
void PField::setGridLocal(double *x_local, double *y_local, double *z_local)
/********************************************************************/
{
	this->x_local_ = x_local;
	this->y_local_ = y_local;
	this->z_local_ = z_local;
}

/********************************************************************/
void PField::setDataOperation( const float *a)
/********************************************************************/
{
	*this = a;
	// long index_local, index_operation;

	// index_operation=0;
	// for (int i = 0; i < this->dims_operation_[0]; i++) {
	// 	for (int j = 0; j < this->dims_operation_[1]; j++) {
	// 		for (int k = 0; k < this->dims_operation_[2]; k++) {
	// 			index_local = this->indexOperationToLocal(i, j, k);
	// 			this->data_local[index_local] = (double)a[index_operation];
	// 			++index_operation;
	// 		}
	// 	}
	// }

	// Set unsynchronized
	this->synchronized_ = false;
}

/********************************************************************/
void PField::setDataOperation( const double *a)
/********************************************************************/
{
	*this = a;
	// long index_local, index_operation;

	// index_operation=0;
	// for (int i = 0; i < this->dims_operation_[0]; i++) {
	// 	for (int j = 0; j < this->dims_operation_[1]; j++) {
	// 		for (int k = 0; k < this->dims_operation_[2]; k++) {
	// 			index_local = this->indexOperationToLocal(i, j, k);
	// 			this->data_local[index_local] = a[index_operation++];
	// 		}
	// 	}
	// }

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
	long index_local, index_operation;

	index_operation=0;
	for (int i = 0; i < this->dims_operation_[0]; i++) {
		for (int j = 0; j < this->dims_operation_[1]; j++) {
			for (int k = 0; k < this->dims_operation_[2]; k++) {
				index_local = this->indexOperationToLocal(i, j, k);
				this->data_local[index_local] = a[index_operation++];
			}
		}
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
PField &PField::operator=( const float *a)
/********************************************************************/
{
	long index_local, index_operation;

	index_operation=0;
	for (int i = 0; i < this->dims_operation_[0]; i++) {
		for (int j = 0; j < this->dims_operation_[1]; j++) {
			for (int k = 0; k < this->dims_operation_[2]; k++) {
				index_local = this->indexOperationToLocal(i, j, k);
				this->data_local[index_local] = (double)a[index_operation++];
			}
		}
	}

	// Set unsynchronized
	this->synchronized_ = false;
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

	long index;

	for (int i = 0; i < this->dims_operation_[0]; i++) {
		for (int j = 0; j < this->dims_operation_[1]; j++) {
			for (int k = 0; k < this->dims_operation_[2]; k++) {
				index = this->indexOperationToLocal(i, j, k);
				this->data_local[index] += a.data_local[index];
			}
		}
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
PField &PField::operator+=( const double *a)
/********************************************************************/
{
	long index_local, index_operation;

	index_operation=0;
	for (int i = 0; i < this->dims_operation_[0]; i++) {
		for (int j = 0; j < this->dims_operation_[1]; j++) {
			for (int k = 0; k < this->dims_operation_[2]; k++) {
				index_local = this->indexOperationToLocal(i, j, k);
				this->data_local[index_local] += a[index_operation++];
			}
		}
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
PField &PField::operator+=( const float *a)
/********************************************************************/
{
	long index_local, index_operation;

	index_operation=0;
	for (int i = 0; i < this->dims_operation_[0]; i++) {
		for (int j = 0; j < this->dims_operation_[1]; j++) {
			for (int k = 0; k < this->dims_operation_[2]; k++) {
				index_local = this->indexOperationToLocal(i, j, k);
				this->data_local[index_local] += (double)a[index_operation++];
			}
		}
	}

	// Set unsynchronized
	this->synchronized_ = false;
	return *this;
}
/********************************************************************/
PField& PField::operator+=(double c)
/********************************************************************/
{

	long index;

	for (int i = 0; i < this->dims_operation_[0]; i++) {
		for (int j = 0; j < this->dims_operation_[1]; j++) {
			for (int k = 0; k < this->dims_operation_[2]; k++) {
				index = this->indexOperationToLocal(i, j, k);
				this->data_local[index] += c;
			}
		}
	}

	// Set unsynchronized
	this->synchronized_ = false;

	return *this;
}

/********************************************************************/
PField& PField::operator-=(const PField &a)
/********************************************************************/
{

	long index;

	for (int i = 0; i < this->dims_operation_[0]; i++) {
		for (int j = 0; j < this->dims_operation_[1]; j++) {
			for (int k = 0; k < this->dims_operation_[2]; k++) {
				index = this->indexOperationToLocal(i, j, k);
				this->data_local[index] -= a.data_local[index];
			}
		}
	}

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
	long index_local, index_operation;

	index_operation=0;
	for (int i = 0; i < this->dims_operation_[0]; i++) {
		for (int j = 0; j < this->dims_operation_[1]; j++) {
			for (int k = 0; k < this->dims_operation_[2]; k++) {
				index_local = this->indexOperationToLocal(i, j, k);
				this->data_local[index_local] -= a[index_operation++];
			}
		}
	}

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
	long index_local, index_operation;

	index_operation=0;
	for (int i = 0; i < this->dims_operation_[0]; i++) {
		for (int j = 0; j < this->dims_operation_[1]; j++) {
			for (int k = 0; k < this->dims_operation_[2]; k++) {
				index_local = this->indexOperationToLocal(i, j, k);
				this->data_local[index_local] -= (double)a[index_operation++];
			}
		}
	}

	// Set unsynchronized
	this->synchronized_ = false;
	return *this;
}

/********************************************************************/
PField& PField::operator-=(double c)
/********************************************************************/
{

	long index;

	for (int i = 0; i < this->dims_operation_[0]; i++) {
		for (int j = 0; j < this->dims_operation_[1]; j++) {
			for (int k = 0; k < this->dims_operation_[2]; k++) {
				index = this->indexOperationToLocal(i, j, k);
				this->data_local[index] -= c;
			}
		}
	}

	// Set unsynchronized
	this->synchronized_ = false;

	return *this;
}

/********************************************************************/
PField& PField::operator*=(const PField &a)
/********************************************************************/
{
	long index;

	for (int i = 0; i < this->dims_operation_[0]; i++) {
		for (int j = 0; j < this->dims_operation_[1]; j++) {
			for (int k = 0; k < this->dims_operation_[2]; k++) {
				index = this->indexOperationToLocal(i, j, k);
				this->data_local[index] *= a.data_local[index];
			}
		}
	}

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
	long index_local, index_operation;

	index_operation=0;
	for (int i = 0; i < this->dims_operation_[0]; i++) {
		for (int j = 0; j < this->dims_operation_[1]; j++) {
			for (int k = 0; k < this->dims_operation_[2]; k++) {
				index_local = this->indexOperationToLocal(i, j, k);
				this->data_local[index_local] *= a[index_operation++];
			}
		}
	}

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
	long index_local, index_operation;

	index_operation=0;
	for (int i = 0; i < this->dims_operation_[0]; i++) {
		for (int j = 0; j < this->dims_operation_[1]; j++) {
			for (int k = 0; k < this->dims_operation_[2]; k++) {
				index_local = this->indexOperationToLocal(i, j, k);
				this->data_local[index_local] *= (double)a[index_operation++];
			}
		}
	}

	// Set unsynchronized
	this->synchronized_ = false;
	return *this;
}

/********************************************************************/
PField& PField::operator*=(double c)
/********************************************************************/
{
	long index;

	for (int i = 0; i < this->dims_operation_[0]; i++) {
		for (int j = 0; j < this->dims_operation_[1]; j++) {
			for (int k = 0; k < this->dims_operation_[2]; k++) {
				index = this->indexOperationToLocal(i, j, k);
				this->data_local[index] *= c;
			}
		}
	}

	// Set unsynchronized
	this->synchronized_ = false;

	return *this;
}

/********************************************************************/
PField& PField::operator/=(const PField &a)
/********************************************************************/
{

	long index;

	for (int i = 0; i < this->dims_operation_[0]; i++) {
		for (int j = 0; j < this->dims_operation_[1]; j++) {
			for (int k = 0; k < this->dims_operation_[2]; k++) {
				index = this->indexOperationToLocal(i, j, k);
				this->data_local[index] /= a.data_local[index];
			}
		}
	}

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
	long index_local, index_operation;

	index_operation=0;
	for (int i = 0; i < this->dims_operation_[0]; i++) {
		for (int j = 0; j < this->dims_operation_[1]; j++) {
			for (int k = 0; k < this->dims_operation_[2]; k++) {
				index_local = this->indexOperationToLocal(i, j, k);
				this->data_local[index_local] /= a[index_operation++];
			}
		}
	}

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
	long index_local, index_operation;

	index_operation=0;
	for (int i = 0; i < this->dims_operation_[0]; i++) {
		for (int j = 0; j < this->dims_operation_[1]; j++) {
			for (int k = 0; k < this->dims_operation_[2]; k++) {
				index_local = this->indexOperationToLocal(i, j, k);
				this->data_local[index_local] /= (double)a[index_operation++];
			}
		}
	}

	// Set unsynchronized
	this->synchronized_ = false;
	return *this;
}

/********************************************************************/
PField& PField::operator/=(double c)
/********************************************************************/
{

	long index;

	for (int i = 0; i < this->dims_operation_[0]; i++) {
		for (int j = 0; j < this->dims_operation_[1]; j++) {
			for (int k = 0; k < this->dims_operation_[2]; k++) {
				index = this->indexOperationToLocal(i, j, k);
				this->data_local[index] /= c;
			}
		}
	}

	// Set unsynchronized
	this->synchronized_ = false;

	return *this;
}

//////////////////////////////////////////////////////////////////////
/// MEMBER FUNCTIONS
//////////////////////////////////////////////////////////////////////

/********************************************************************/
PField &PField::add(PField &a, PField &b)
/********************************************************************/
{

#ifdef BOUNDS_CHECK
	long N = this->getSizeOperation();
	// First check that the fields are the same size
	if ( N != a.getSizeOperation() || N != b.getSizeOperation() ) {
		cout << "Mismatch in field sizes" << endl;
		exit (EXIT_FAILURE);
	}
#endif

	long index;
	for (int i = 0; i < this->dims_operation_[0]; i++) {
		for (int j = 0; j < this->dims_operation_[1]; j++) {
			for (int k = 0; k < this->dims_operation_[2]; k++) {
				index = this->indexOperationToLocal(i, j, k);
				this->data_local[index] = a.data_local[index]
					+ b.data_local[index];
			}
		}
	}

	// Set unsynchronized
	this->synchronized_ = false;
	return *this;
}

/********************************************************************/
PField &PField::sub(PField &a, PField &b)
/********************************************************************/
{

#ifdef BOUNDS_CHECK
	long N = this->getSizeOperation();
	// First check that the fields are the same size
	if ( N != a.getSizeOperation() || N != b.getSizeOperation() ) {
		cout << "Mismatch in field sizes" << endl;
		exit (EXIT_FAILURE);
	}
#endif

	long index;
	for (int i = 0; i < this->dims_operation_[0]; i++) {
		for (int j = 0; j < this->dims_operation_[1]; j++) {
			for (int k = 0; k < this->dims_operation_[2]; k++) {
				index = this->indexOperationToLocal(i, j, k);
				this->data_local[index] = a.data_local[index]
					- b.data_local[index];
			}
		}
	}

	// Set unsynchronized
	this->synchronized_ = false;
	return *this;
}

/********************************************************************/
PField &PField::mul(PField &a, PField &b)
/********************************************************************/
{

#ifdef BOUNDS_CHECK
	long N = this->getSizeOperation();
	// First check that the fields are the same size
	if ( N != a.getSizeOperation() || N != b.getSizeOperation() ) {
		cout << "Mismatch in field sizes" << endl;
		exit (EXIT_FAILURE);
	}
#endif

	long index;
	for (int i = 0; i < this->dims_operation_[0]; i++) {
		for (int j = 0; j < this->dims_operation_[1]; j++) {
			for (int k = 0; k < this->dims_operation_[2]; k++) {
				index = this->indexOperationToLocal(i, j, k);
				this->data_local[index] = a.data_local[index]
					* b.data_local[index];
			}
		}
	}

	// Set unsynchronized
	this->synchronized_ = false;
	return *this;
}

/********************************************************************/
PField &PField::div(PField &a, PField &b)
/********************************************************************/
{

#ifdef BOUNDS_CHECK
	long N = this->getSizeOperation();
	// First check that the fields are the same size
	if ( N != a.getSizeOperation() || N != b.getSizeOperation() ) {
		cout << "Mismatch in field sizes" << endl;
		exit (EXIT_FAILURE);
	}
#endif

	long index;
	for (int i = 0; i < this->dims_operation_[0]; i++) {
		for (int j = 0; j < this->dims_operation_[1]; j++) {
			for (int k = 0; k < this->dims_operation_[2]; k++) {
				index = this->indexOperationToLocal(i, j, k);
				this->data_local[index] = a.data_local[index]
					/ b.data_local[index];
			}
		}
	}

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

	void (FiniteDiff::*dd)(int, int, double *, int, double *);
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

	void (FiniteDiff::*dd)(int, int, double *, int, double *);
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

	void (FiniteDiff::*dd)(int, int, double *, int, double *);
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

	void (FiniteDiff::*dd)(int, int, double *, int, double *);
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

	void (FiniteDiff::*dd)(int, int, double *, int, double *);
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

	void (FiniteDiff::*dd)(int, int, double *, int, double *);
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

	void (FiniteDiff::*dd)(int, int, double *, int, double *);
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

	void (FiniteDiff::*dd)(int, int, double *, int, double *);
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

	void (FiniteDiff::*dd)(int, int, double *, int, double *);
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
void PField::dndxn(void (FiniteDiff::*dd)(int, int, double *, int, double *), PField &a)
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
void PField::dndyn(void (FiniteDiff::*dd)(int, int, double *, int, double *),
		  PField &a)
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
void PField::dndzn(void (FiniteDiff::*dd)(int, int, double *, int, double *),
		  PField &a)
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
int *PField::computeMPITopologyDims(int nproc, int mpi_decomp_ndims)
/********************************************************************/
{

	// initialize the MPI topology dimensions to 1
	int *mpi_topology_dims = new int[PFIELD_NDIMS];
	std::fill_n(mpi_topology_dims, PFIELD_NDIMS, 1);

	if (mpi_decomp_ndims == 1) {

		mpi_topology_dims[0] = nproc;

	} else if (mpi_decomp_ndims == 2) {

		int M = std::min( (int) ceil( sqrt(nproc * (double) this->dims_[0] / this->dims_[1])), nproc);
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
bool PField::hasRind(int dim, int location)
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

	const long prev_buffer_size = this->getSizeRind(dim, -1);
	const long next_buffer_size = this->getSizeRind(dim, 1);
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
double *PField::createRindBuffer(int dim, int location)
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
void PField::packRindBuffer(int dim, int location, double *rind_buffer)
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

	long index = 0;
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
void PField::unpackRindBuffer(int dim, int location, double *rind_buffer)
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

	// Reset values along dim direction
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

	long index = 0;
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

