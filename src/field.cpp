#include <iostream>
#include <cstdlib>
#include <string.h>
#include "field.hpp"

using namespace std;

namespace pturb_fields {

  ////////////////////////////////////////////////////////////////////////////////
  /// CONSTRUCTORS
  //////////////////////////////////////////////////////////////////////////////// 
  
    // Field class member functions
  /********************************************************************/
  Field::Field( MPI_Comm comm, FieldDecomp_t field_decomp, const int *dims);
  /********************************************************************/
  {
    cout << "Initializing Field" << endl;
    FieldInit( comm, decomp, ndim, dims );
  }

  /********************************************************************/
  Field::Field( const Field &g ) 
  /********************************************************************/
  {
    cout << "Copy Constructor Initializing Field" << endl;

    FieldInit( g.getMpiComm(), g.getFieldDecomp(), g.getDims() );

    // Copy the field data from g to this->data
    memcpy(this->data_,g.data_,sizeof(*this->data_)*this->getSize());

  }

  ////////////////////////////////////////////////////////////////////////////////
  /// DECONSTRUCTORS
  ////////////////////////////////////////////////////////////////////////////////

  /********************************************************************/
  Field::~Field(void) 
  /********************************************************************/
  {
    cout << "Field deconstructor" << endl;
    delete[] this->data_;
    delete this->finite_diff_;
    delete this-> mpi_topology_;
  }

  ////////////////////////////////////////////////////////////////////////////////
  /// INITIALIZERS
  //////////////////////////////////////////////////////////////////////////////// 

  // Field class initializer
  // Sets:
  //   ndims
  // Initilizes:
  //   dims
  //   data
  /********************************************************************/
  void Field::FieldInit( int *dims, int *periodic, FieldDecomp_t field_decomp, int operator_order ) 
  /********************************************************************/
  {

  
    // Set the global dimensions
    memcpy( this->dims_, dims, sizeof(*this->dims_) * FIELD_NDIMS );
    

    int mpi_coord_ndims;
    int *mpi_coord_dims;

    // Determine the MPI coordinate dimensions
    if( field_decomp == FIELD_DECOMP_SLAB ) {
      mpi_coord_ndims = 1;
      


    // Create the MPI topology struct
    this->mpi_topology_ = MpiTopologyNew( MPI_COMM_WORLD, mpi_coord_ndims, mpi_coord_dims, periodic );

    // Get the processor local dimensions from the domain
    // decomposition; may have struct here for additional data
    int *dims_local = domainDecomp( _nproc, _rank, _decomp, _dims );



    // Initialize the data field
    long N = getSizeLocal();
    this->data_local = new double[N];

    // Initialize pointers to null
    this->finite_diff_ = NULL;

  }

  // Field class initializer; expects ndims to be set
  // Initilizes:
  //   dims
  //   data
  /********************************************************************/
  void Field::FieldInit( int *dims ) 
  /********************************************************************/
  {
    FieldInit( this->ndims, dims );
  }

  // Initialize the finite difference derivatives class
  /********************************************************************/
  void Field::derivFDInit( int order )
  /********************************************************************/
  {

    if( this->finite_diff_ != NULL ) {
      cout << "finite difference class instance fd already allocated" << endl;
      exit (EXIT_FAILURE);
    }

    // First initialize the fd class instance
    this->finite_diff_ = new FiniteDiff( this->dims_local_[0], this->x_local_, 
					 this->dims_local_[1], this->y_local_, 
					 this->dims_local_[2], this->z_local_, 
					 order );
    
  }

  /********************************************************************/
  int *Field::domainDecomp( int nproc, int rank, FieldDecomp_t field_decomp, const int *dims )
  /********************************************************************/
  {

    int dims_local[FIELD_NDIMS];

    // Now need to decompose 
    if( field_decomp == FIELD_DECOMP_SLAB ) {
      // Cut along the first dimension
      dims_loc[0] = dims[0] / nproc;
      for (int n=1; n<FIELD_NDIMS; n++ ) dims_local[n] = dims[n];
    } else {
      cout << "Only FieldDecomp_t FIELD_DECOMP_SLAB currently supported" << endl;
      exit (EXIT_FAILURE);
    }    
    
    return dims_local;

  }

  ////////////////////////////////////////////////////////////////////// 
  /// GETTERS/SETTERS
  //////////////////////////////////////////////////////////////////////

  /********************************************************************/
  long Field::getSize()
  /********************************************************************/
  {
    long N=1;
    for (int n=0; n<this->FIELD_NDIMS; n++ ) {
      N *= this->dims_[n];
    }
    return N;
  }

  /********************************************************************/
  MPI_Comm Field::getMpiComm()
  /********************************************************************/
  {
    return this->mpi_topology_->comm;
  }

  /********************************************************************/
  FieldDecomp_t Field::getFieldDecomp()
  /********************************************************************/
  {
    return this->field_decomp_;
  }

  /********************************************************************/
  int *Field::getDims()
  /********************************************************************/
  {
    return this->dims_;
  }

  // Array index for 3D indexing
  /********************************************************************/
  long Field::index( int i, int j, int k )
  /********************************************************************/
  {
    return ((long)this->dims_[1]*(long)i + (long)j)*(long)this->dims_[2] + (long)k;
  }

  // Array index for 3D indexing
  /********************************************************************/
  long Field::indexLocal( int i, int j, int k )
  /********************************************************************/
  {
    return ((long)this->dims_local_[1]*(long)i + (long)j)*(long)this->dims_local_[2] + (long)k;
  }

  // Array index for 3D indexing
  /********************************************************************/
  long Field::indexOperation( int i, int j, int k )
  /********************************************************************/
  {
    return ((long)this->dims_operation_[1]*(long)i + (long)j)*(long)this->dims_operation_[2] + (long)k;
  }

  // Array index for 3D indexing
  /********************************************************************/
  long Field::indexOperationToLocal( int i, int j, int k )
  /********************************************************************/
  {
    static const long ii = (long)i + (long)this->offset_operation_[0];
    static const long jj = (long)j + (long)this->offset_operation_[1];
    static const long kk = (long)k + (long)this->offset_operation_[2];

    return ((long)this->dims_local_[1]*ii + jj) * (long)this->dims_local_[2] + kk;
  }

  // Set the grid pointers for field
  /********************************************************************/
  void Field::assignGridLocal( double *x_local, double *y_local, double *z_local )
  /********************************************************************/
  {  
    this->x_local_ = x_local;
    this->y_local_ = y_local;
    this->z_local_ = z_local;
  }

  //////////////////////////////////////////////////////////////////////
  /// OPERATORS
  ////////////////////////////////////////////////////////////////////// 

  /********************************************************************/
  Field& Field::operator=( const Field &a )
  /********************************************************************/
  {
    // Make sure they don't point to the same address
    if (this != &a) {
      // Copy the field data
      memcpy( this->data_local, a.data_local, sizeof(*this->data_local) * this->getSizeLocal() );
    }	
    return *this;
  }


  /********************************************************************/
  Field& Field::operator+=( const Field &a )
  /********************************************************************/
  {

    long index;

    for( int i=0; i<this->dims_operation_[0]; i++ ) {
      for (int j=0; j<this->dims_operation_[1]; j++ ) {
	for ( int k=0; k<this->dims_operation_[2]; k++ ) {
	  index = this->indexOperationToLocal(i,j,k);
	  this->data_local[index] += a.data_local[index];
	}
      }
    }

    return *this;
  }

  /********************************************************************/
  Field& Field::operator-=( const Field &a )
  /********************************************************************/
  {

    long index;

    for( int i=0; i<this->dims_operation_[0]; i++ ) {
      for (int j=0; j<this->dims_operation_[1]; j++ ) {
	for ( int k=0; k<this->dims_operation_[2]; k++ ) {
	  index = this->indexOperationToLocal(i,j,k);
	  this->data_local[index] -= a.data_local[index];
	}
      }
    }
    return *this;
  }

  /********************************************************************/
  Field& Field::operator*=( const Field &a )
  /********************************************************************/
  {
    long index;

    for( int i=0; i<this->dims_operation_[0]; i++ ) {
      for (int j=0; j<this->dims_operation_[1]; j++ ) {
	for ( int k=0; k<this->dims_operation_[2]; k++ ) {
	  index = this->indexOperationToLocal(i,j,k);
	  this->data_local[index] *= a.data_local[index];
	}
      }
    }
	
    return *this;
  }

  /********************************************************************/
  Field& Field::operator/=( const Field &a )
  /********************************************************************/
  {

    long index;

    for( int i=0; i<this->dims_operation_[0]; i++ ) {
      for (int j=0; j<this->dims_operation_[1]; j++ ) {
	for ( int k=0; k<this->dims_operation_[2]; k++ ) {
	  index = this->indexOperationToLocal(i,j,k);
	  this->data_local[index] /= a.data_local[index];
	}
      }
    }
    return *this;
  }


  //////////////////////////////////////////////////////////////////////
  /// MEMBER FUNCTIONS
  //////////////////////////////////////////////////////////////////////


  /********************************************************************/
  void Field::add( Field &a, Field &b )
  /********************************************************************/
  {
 
#ifdef BOUNDS_CHECK
    long N = this->getSizeOperation();
    // First check that the fields are the same size
    if ( N != a.getSizeOperation() || N != b.getSizOperation() ) {
      cout << "Mismatch in field sizes" << endl;
      exit (EXIT_FAILURE);
    }   
#endif

    long index;
    for( int i=0; i<this->dims_operation_[0]; i++ ) {
      for (int j=0; j<this->dims_operation_[1]; j++ ) {
	for ( int k=0; k<this->dims_operation_[2]; k++ ) {
	  index = this->indexOperationToLocal(i,j,k);
	  this->data_local[index] = a.data_local[index] + b.data_local[index];
	}
      }
    }
  }

  /********************************************************************/
  void Field::sub( Field &a, Field &b )
  /********************************************************************/
  {
 
#ifdef BOUNDS_CHECK
    long N = this->getSizeOperation();
    // First check that the fields are the same size
    if ( N != a.getSizeOperation() || N != b.getSizOperation() ) {
      cout << "Mismatch in field sizes" << endl;
      exit (EXIT_FAILURE);
    }   
#endif

    long index;
    for( int i=0; i<this->dims_operation_[0]; i++ ) {
      for (int j=0; j<this->dims_operation_[1]; j++ ) {
	for ( int k=0; k<this->dims_operation_[2]; k++ ) {
	  index = this->indexOperationToLocal(i,j,k);
	  this->data_local[index] = a.data_local[index] - b.data_local[index];
	}
      }
    }
  }

  /********************************************************************/
  void Field::mul( Field &a, Field &b )
  /********************************************************************/
  {
 
#ifdef BOUNDS_CHECK
    long N = this->getSizeOperation();
    // First check that the fields are the same size
    if ( N != a.getSizeOperation() || N != b.getSizOperation() ) {
      cout << "Mismatch in field sizes" << endl;
      exit (EXIT_FAILURE);
    }   
#endif

    long index;
    for( int i=0; i<this->dims_operation_[0]; i++ ) {
      for (int j=0; j<this->dims_operation_[1]; j++ ) {
	for ( int k=0; k<this->dims_operation_[2]; k++ ) {
	  index = this->indexOperationToLocal(i,j,k);
	  this->data_local[index] = a.data_local[index] * b.data_local[index];
	}
      }
    }
  }

  /********************************************************************/
  void Field::div( Field &a, Field &b )
  /********************************************************************/
  {
 
#ifdef BOUNDS_CHECK
    long N = this->getSizeOperation();
    // First check that the fields are the same size
    if ( N != a.getSizeOperation() || N != b.getSizOperation() ) {
      cout << "Mismatch in field sizes" << endl;
      exit (EXIT_FAILURE);
    }   
#endif

    long index;
    for( int i=0; i<this->dims_operation_[0]; i++ ) {
      for (int j=0; j<this->dims_operation_[1]; j++ ) {
	for ( int k=0; k<this->dims_operation_[2]; k++ ) {
	  index = this->indexOperationToLocal(i,j,k);
	  this->data_local[index] = a.data_local[index] / b.data_local[index];
	}
      }
    }
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
  void Field::ddx( Field &a ) 
  /********************************************************************/
  {
    void (FiniteDiff::*dd)(int, double *, double *);
    dd = &FiniteDiff::ddx;
    this->dndxn( dd, a );
  }

  /********************************************************************/
  void Field::ddy( Field &a ) 
  /********************************************************************/
  {
    void (FiniteDiff::*dd)(int, double *, double *);
    dd = &FiniteDiff::ddy;
    this->dndyn( dd, a );
  }

  /********************************************************************/
  void Field::ddz( Field &a ) 
  /********************************************************************/
  {
    void (FiniteDiff::*dd)(int, double *, double *);
    dd = &FiniteDiff::ddz;
    this->dndzn( dd, a );
  }

  /********************************************************************/
  void Field::d2dx2( Field &a ) 
  /********************************************************************/
  {
    void (FiniteDiff::*dd)(int, double *, double *);
    dd = &FiniteDiff::d2dx2;
    this->dndxn( dd, a );
  }
  /********************************************************************/
  void Field::d2dy2( Field &a ) 
  /********************************************************************/
  {
    void (FiniteDiff::*dd)(int, double *, double *);
    dd = &FiniteDiff::d2dy2;
    this->dndyn( dd, a );
  }

  /********************************************************************/
  void Field::d2dz2( Field &a ) 
  /********************************************************************/
  {
    void (FiniteDiff::*dd)(int, double *, double *);
    dd = &FiniteDiff::d2dz2;
    this->dndzn( dd, a );
  }

  /********************************************************************/
  void Field::d2dxy( Field &a ) 
  /********************************************************************/
  {
    void (FiniteDiff::*dd)(int, double *, double *);
    dd = &FiniteDiff::ddx;
    this->dndxn( dd, a );
    dd = &FiniteDiff::ddy;
    this->dndyn( dd, *this );
  } 
  /********************************************************************/
  void Field::d2dxz( Field &a ) 
  /********************************************************************/
  {
    void (FiniteDiff::*dd)(int, double *, double *);
    dd = &FiniteDiff::ddx;
    this->dndxn( dd, a );
    dd = &FiniteDiff::ddz;
    this->dndzn( dd, *this );
  }
  /********************************************************************/
  void Field::d2dyz( Field &a ) 
  /********************************************************************/
  {
    void (FiniteDiff::*dd)(int, double *, double *);
    dd = &FiniteDiff::ddy;
    this->dndyn( dd, a );
    dd = &FiniteDiff::ddz;
    this->dndzn( dd, *this );
  }    

  ////////////////////////////////////////////////////////////////////// 
  /// DERIVATIVES (PRIVATE)
  //////////////////////////////////////////////////////////////////////

  /********************************************************************/
  void Field::dndxn( void (FiniteDiff::*dd)( int, double *, double *), Field &a ) 
  /********************************************************************/
  {

    // Create buffer to store x-data
    const int nx = this->dims_operation_[0];
    const int ny = this->dims_operation_[1];
    const int nz = this->dims_operation_[2];
    
#ifdef BOUNDS_CHECK
    // First check that the fields are the same size
    if ( nx != a.dims_opteration_[0] || ny != a.dims_operation_[1] || nz != a.dims_operation_[2] ) {
      cout << "Mismatch in field sizes" << endl;
      exit (EXIT_FAILURE);
    }   
#endif

    double *ax = new double[nx];
    double *dax = new double[nx];
    
    for (int j=0; j<ny; j++ ) {
      for (int k=0; k<nz; k++ ) {
	// Pack buffer
	for ( int i=0; i<nx; i++ ) ax[i] = a.data_local[a.indexOperationToLocal(i,j,k)];
	// Compute the derivative using the FiniteDiff class
	(this->finite_diff_->*dd)( this->offset_operation_[0], nx, ax, dax );
	// Unpack buffer
	for ( int i=0; i<nx; i++ ) this->data_local[this->indexOperationToLocal(i,j,k)] = dax[i];
      }
    }

    delete [] ax;
    delete [] dax;
    
  }  

  /********************************************************************/
  void Field::dndyn( void (FiniteDiff::*dd)( int, double *, double *), Field &a ) 
  /********************************************************************/
  {

    // Create buffer to store x-data
    const int nx = this->dims_operation_[0];
    const int ny = this->dims_operation_[1];
    const int nz = this->dims_operation_[2];
    
#ifdef BOUNDS_CHECK
    // First check that the fields are the same size
    if ( nx != a.dims_opteration_[0] || ny != a.dims_operation_[1] || nz != a.dims_operation_[2] ) {
      cout << "Mismatch in field sizes" << endl;
      exit (EXIT_FAILURE);
    }   
#endif

    double *ay = new double[ny];
    double *day = new double[ny];
    
    for ( int i=0; i<nx; i++ ) {
      for (int k=0; k<nz; k++ ) {
	// Pack buffer
	for (int j=0; j<ny; j++ ) ay[j] = a.data_local[a.indexOperationToLocal(i,j,k)];
	// Compute the derivative using the FiniteDiff class
	(this->finite_diff_->*dd)( this->offset_operation_[1], ny, ay, day );
	// Unpack buffer
	for ( int j=0; j<ny; j++ ) this->data_[this->indexOperationToLocal(i,j,k)] = day[j];
      }
    }

    delete [] ay;
    delete [] day;
    
  }

  /********************************************************************/
  void Field::dndzn( void (FiniteDiff::*dd)( int, double *, double *), Field &a ) 
  /********************************************************************/
  {

    // Create buffer to store x-data
    const int nx = this->dims_operation_[0];
    const int ny = this->dims_operation_[1];
    const int nz = this->dims_operation_[2];
    
#ifdef BOUNDS_CHECK
    // First check that the fields are the same size
    if ( nx != a.dims_opteration_[0] || ny != a.dims_operation_[1] || nz != a.dims_operation_[2] ) {
      cout << "Mismatch in field sizes" << endl;
      exit (EXIT_FAILURE);
    }   
#endif

    double *az = new double[nz];
    double *daz = new double[nz];
    
    for ( int i=0; i<nx; i++ ) {
      for (int j=0; j<ny; j++ ) {
	// Pack buffer
	for (int k=0; k<nz; k++ ) az[k] = a.data_local[a.indexOperationToLocal(i,j,k)];
	// Compute the derivative using the FiniteDiff class
	(this->finite_diff_->*dd)( this->offset_operation_[2], nz, az, daz );
	// Unpack buffer
	for (int k=0; k<nz; k++ ) this->data_local[this->indexOperationToLocal(i,j,k)] = daz[k];
      }
    }

    delete [] az;
    delete [] daz;
    
  }


  //////////////////////////////////////////////////////////////////////
  /// NON-CLASS MEMBER FUNCTIONS
  //////////////////////////////////////////////////////////////////////

}


