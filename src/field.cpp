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
  Field::Field( const int *dims, FieldDecomp_t field_decomp, const int *periodic, int operator_order ) 
  /********************************************************************/
  {
    cout << "Initializing Field" << endl;
    FieldInit( dims, field_decomp, periodic, operator_order );
  }

  /********************************************************************/
  Field::Field( const Field &g ) 
  /********************************************************************/
  {
    cout << "Copy Constructor Initializing Field" << endl;

    FieldInit( g.getDims(), g.getFieldDecomp(), g.getMpiTopology()->periodic, g.getOperatorOrder() );

    // Copy the field data from g to this->data
    memcpy(this->data_local,g.data_local,sizeof(*this->data_local)*this->getSizeLocal());

  }

  ////////////////////////////////////////////////////////////////////////////////
  /// DECONSTRUCTORS
  ////////////////////////////////////////////////////////////////////////////////

  /********************************************************************/
  Field::~Field(void) 
  /********************************************************************/
  {
    cout << "Field deconstructor" << endl;
    delete[] this->data_local;
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
  void Field::FieldInit( int *dims, FieldDecomp_t field_decomp, int *periodic, int operator_order ) 
  /********************************************************************/
  {

  
    // Set the global dimensions
    memcpy( this->dims_, dims, sizeof(*this->dims_) * FIELD_NDIMS );
    this->field_decomp_ = field_decomp;
    this->operator_order_ = operator_order;

    int mpi_coord_ndims;
    int *mpi_coord_dims;

    // Determine the MPI coordinate dimensions
    if( field_decomp == FIELD_DECOMP_SLAB ) {
      mpi_coord_ndims = 1;
    } else if ( field_decomp == FIELD_DECOMP_PENCIL ) {
      mpi_coord_ndims = 2;
    } else {
      cout << "only slab or pencil decompositions currently supported" << endl;
      exit(EXIT_FAILURE);
    }
    // Compute the coordinate dimensions for the given Field decomposition
    mpi_coord_dims = computeMpiTopologyDims( mpi_coord_ndims );

    // Create the MPI topology struct; this contains all the comm, rank, nproc, neighbors, etc.
    this->mpi_topology_ = MpiTopologyNew( MPI_COMM_WORLD, mpi_coord_ndims, mpi_coord_dims, periodic );

    int *zero3=new int(0)[FIELD_NDIMS]; // Zero vector of dimension FIELD_NDIMS

    // Rind size to support finite differencing at interprocessor boundaries
    int rind_size = operator_order / 2 + operator_order % 2; // Second part adds 1 if odd operator order

    // Initialize the local/operation dimensions and offset
    memcpy( this->dims_local_, this->dims_, sizeof( *this->dims_ ) * FIELD_NDIMS );
    memcpy( this->offset_local_, zero3, sizeof( *this->offset_local_ ) * FIELD_NDIMS );

    memcpy( this->dims_operation_, this->dims_, sizeof( *this->dims_ ) * FIELD_NDIMS );
    memcpy( this->offset_operation_, zero3, sizeof( *this->offset_local_ ) * FIELD_NDIMS );

    delete [] zero3;

    // For each dimension in the MPI topology; determin what the dimensions and
    // offsets are of the local and operation domains
    for( int n=0; n<this->mpi_topology_->ndims; n++ ) {

      int p = this->dims_[n] % this->mpi_topology_->dims[n]; // Remainder of points (surplus)

      int d = this->dims_[n] / this->mpi_topology_->dims[n];
      int s = ( this->mpi_topology_->coords[n] < p ) ? 1 : 0; // Take one of the surplus values if in low enough rank

      // Check that the rind size is not larger than the size of the smallest dimension
      if( d < rind_size ) {
	cout << "Domain decomposition too fine to support derivatives of order " << operator_order << endl;
	exit(EXIT_FAILURE);
      }

      this->dims_local_[n] = d + s; // Local domain size without rind data
      // Offset of the local domain with respect to the global domain
      this->offset_local_[n] = this->mpi_topology_->coords[n] * d + std::min( this->mpi_topology_->coords[n], p );

      // The operation domain is the local domain without rind data
      this->dims_operation_[n] = this->dims_local_[n];

      // Now determine the full local dimensions with the rind
      if( this->mpi_topology_->neighbors_next[n] >= 0 ) this->dims_local_[n] += rind_size; // If we have a neighbor then add the rind points
      if( this->mpi_topology_->neighbors_prev[n] >= 0 ) {
	this->dims_local_[n] += rind_size;
	this->offset_local_[n] -= rind_size; // Move the local offset position back to compensate for the rind points
	this->offset_operation_[n] = rind_size; // Set the operation offset distance relative to the local domain 
      }

    }

    // Now that we have set dims_local_ we may get the total number of elements
    this->data_local = new double[this->getSizeLocal()];

    // Initialize pointers to null
    this->finite_diff_ = NULL;

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
  int *Field::computeMpiTopologyDims( mpi_coord_ndims ) 
  /********************************************************************/
  {

    int nproc;

    // Get the total number of procs
    MPI_Comm_size( MPI_COMM_WORLD, &nproc );


    int *mpi_topology_dims = new int[mpi_coord_ndims];
    if( mpi_coord_ndims == 1 ) {
      mpi_topology_dims[0] = nproc;
    } else if( mpi_coord_ndims == 2 ) {

	int M = std::min( (int)ceil( sqrt(nproc*(double)this->dims_[0]/this->dims_[1] )), nproc );
	int N = nproc / M;
	
	while(1) {
	  if( N*M == nproc ) {
	    break;
	  } else if( N*M < nproc ) {
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
  MpiTopology_t *Field::getMpiTopology()
  /********************************************************************/
  {
    return this->mpi_topology_;
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


