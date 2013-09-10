#include <iostream>
#include <cstdlib>
#include <string.h>
#include "field.hpp"

using namespace std;

namespace pturb_fields {

    // Field class member functions
  /********************************************************************/
  Field::Field( int ndim, int *dims ) 
  /********************************************************************/
  {
    cout << "Initializing Field" << endl;
    FieldInit( ndim, dims );
  }

  /********************************************************************/
  Field::Field( const Field &g ) 
  /********************************************************************/
  {
    cout << "Copy Constructor Initializing Field" << endl;

    FieldInit( g.ndim_, g.dims_ );
    // Copy the field data from g to this->data
    memcpy(this->data_,g.data_,sizeof(*this->data_)*this->getSize());

  }

  /********************************************************************/
  Field::~Field(void) 
  /********************************************************************/
  {
    cout << "Field deconstructor" << endl;
    delete[] this->dims_;
    delete[] this->data_;
  }

  // Field class initializer; expects ndim to already be set.
  // Initilizes:
  //   dims
  //   data
  /********************************************************************/
  void Field::FieldInit( int *dims ) 
  /********************************************************************/
  {
    this->dims_ = new int [this->ndim_];    
    // Copy the dims vector
    memcpy( this->dims_, dims, sizeof(*this->dims_) * this->ndim_ );

    // Initialize the data field
    long N = getSize();
    this->data_ = new double[N];

    // Initialize pointers to null
    this->fd_ = NULL;
  }

  // Field class initializer
  // Sets:
  //   ndim
  // Initilizes:
  //   dims
  //   data
  /********************************************************************/
  void Field::FieldInit( int ndim, int *dims ) 
  /********************************************************************/
  {
    this->ndim_ = ndim;
    FieldInit( dims );
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
      memcpy( this->data_, a.data_, sizeof(*this->data_) * this->getSize() );
    }	
    return *this;
  }


  /********************************************************************/
  Field& Field::operator+=( const Field &a )
  /********************************************************************/
  {
    const long N=this->getSize();
    long indx=0;
    while( indx != N ) {
      this->data_[indx] += a.data_[indx];
      ++indx;
    }	
    return *this;
  }

  /********************************************************************/
  Field& Field::operator-=( const Field &a )
  /********************************************************************/
  {
    const long N=this->getSize();
    long indx=0;
    while( indx != N ) {
      this->data_[indx] -= a.data_[indx];
      ++indx;
    }	
    return *this;
  }

  /********************************************************************/
  Field& Field::operator*=( const Field &a )
  /********************************************************************/
  {
    const long N=this->getSize();
    long indx=0;
    while( indx != N ) {
      this->data_[indx] *= a.data_[indx];
      ++indx;
    }	
    return *this;
  }

  /********************************************************************/
  Field& Field::operator/=( const Field &a )
  /********************************************************************/
  {
    const long N=this->getSize();
    long indx=0;
    while( indx != N ) {
      this->data_[indx] /= a.data_[indx];
      ++indx;
    }	
    return *this;
  }


  //////////////////////////////////////////////////////////////////////
  /// MEMBER FUNCTIONS
  //////////////////////////////////////////////////////////////////////

  /********************************************************************/
  long Field::getSize()
  /********************************************************************/
  {
    long N=1;
    for (int n=0; n<this->ndim_; n++ ) {
      N *= this->dims_[n];
    }
    return N;
  }

  /********************************************************************/
  void Field::add( Field &a, Field &b )
  /********************************************************************/
  {
    long N = this->getSize();

#ifdef BOUNDS_CHECK
    // First check that the fields are the same size
    if ( N != a.getSize() || N != b.getSize() ) {
      cout << "Mismatch in field sizes" << endl;
      exit (EXIT_FAILURE);
    }   
#endif

    // Perform the multplication
    long indx=0;
    while( indx != N ) {
      this->data_[indx] = a.data_[indx] + b.data_[indx];
      ++indx;
    }
  }

  /********************************************************************/
  void Field::sub( Field &a, Field &b )
  /********************************************************************/
  {
    long N = this->getSize();

#ifdef BOUNDS_CHECK
    // First check that the fields are the same size
    if ( N != a.getSize() || N != b.getSize() ) {
      cout << "Mismatch in field sizes" << endl;
      exit (EXIT_FAILURE);
    }   
#endif
    // Perform the multplication
    long indx=0;
    while( indx != N ) {
      this->data_[indx] = a.data_[indx] - b.data_[indx];
      ++indx;
    }
  }

  // Multiply two fields the resulting fields as this = a*b
  /********************************************************************/
  void Field::mul( Field &a, Field &b ) 
  /********************************************************************/
  {
    long N = this->getSize();

#ifdef BOUNDS_CHECK
    // First check that the fields are the same size
    if ( N != a.getSize() || N != b.getSize() ) {
      cout << "Mismatch in field sizes" << endl;
      exit (EXIT_FAILURE);
    }   
#endif
    // Perform the multplication
    long indx=0;
    while( indx != N ) {
      this->data_[indx] = a.data_[indx] * b.data_[indx];
      ++indx;
    }
  }

  /********************************************************************/
  void Field::div( Field &a, Field &b ) 
  /********************************************************************/
  {
    long N = this->getSize();

#ifdef BOUNDS_CHECK
    // First check that the fields are the same size
    if ( N != a.getSize() || N != b.getSize() ) {
      cout << "Mismatch in field sizes" << endl;
      exit (EXIT_FAILURE);
    }   
#endif
    // Perform the multplication
    long indx=0;
    while( indx != N ) {
      this->data_[indx] = a.data_[indx] / b.data_[indx];
      ++indx;
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

  /********************************************************************/
  void Field::dndxn( void (FiniteDiff::*dd)( int, double *, double *), Field &a ) 
  /********************************************************************/
  {

    // Create buffer to store x-data
    const int nx = this->dims_[0];
    const int ny = this->dims_[1];
    const int nz = this->dims_[2];
    
#ifdef BOUNDS_CHECK
    // First check that the fields are the same size
    if ( nx != a.dims_[0] || ny != a.dims_[1] || nz != a.dims_[2] ) {
      cout << "Mismatch in field sizes" << endl;
      exit (EXIT_FAILURE);
    }   
#endif

    double *ax = new double[nx];
    double *dax = new double[nx];
    
    for (int j=0; j<ny; j++ ) {
      for (int k=0; k<nz; k++ ) {
	// Pack buffer
	for ( int i=0; i<nx; i++ ) ax[i] = a.data_[a.index(i,j,k)];
	// Compute the derivative using the FiniteDiff class
	(this->fd_->*dd)( nx, ax, dax );
	// Unpack buffer
	for ( int i=0; i<nx; i++ ) this->data_[this->index(i,j,k)] = dax[i];
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
    const int nx = this->dims_[0];
    const int ny = this->dims_[1];
    const int nz = this->dims_[2];
    
#ifdef BOUNDS_CHECK
    // First check that the fields are the same size
    if ( nx != a.dims_[0] || ny != a.dims_[1] || nz != a.dims_[2] ) {
      cout << "Mismatch in field sizes" << endl;
      exit (EXIT_FAILURE);
    }   
#endif

    double *ay = new double[ny];
    double *day = new double[ny];
    
    for ( int i=0; i<nx; i++ ) {
      for (int k=0; k<nz; k++ ) {
	// Pack buffer
	for (int j=0; j<ny; j++ ) ay[j] = a.data_[a.index(i,j,k)];
	// Compute the derivative using the FiniteDiff class
	(this->fd_->*dd)( ny, ay, day );
	// Unpack buffer
	for ( int j=0; j<ny; j++ ) this->data_[this->index(i,j,k)] = day[j];
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
    const int nx = this->dims_[0];
    const int ny = this->dims_[1];
    const int nz = this->dims_[2];
    
#ifdef BOUNDS_CHECK
    // First check that the fields are the same size
    if ( nx != a.dims_[0] || ny != a.dims_[1] || nz != a.dims_[2] ) {
      cout << "Mismatch in field sizes" << endl;
      exit (EXIT_FAILURE);
    }   
#endif

    double *az = new double[nz];
    double *daz = new double[nz];
    
    for ( int i=0; i<nx; i++ ) {
      for (int j=0; j<ny; j++ ) {
	// Pack buffer
	for (int k=0; k<nz; k++ ) az[k] = a.data_[a.index(i,j,k)];
	// Compute the derivative using the FiniteDiff class
	(this->fd_->*dd)( nz, az, daz );
	// Unpack buffer
	for (int k=0; k<nz; k++ ) this->data_[this->index(i,j,k)] = daz[k];
      }
    }

    delete [] az;
    delete [] daz;
    
  }


  // Array index for 3D indexing
  /********************************************************************/
  long Field::index( int i, int j, int k )
  /********************************************************************/
  {
    return ((long)this->dims_[1]*(long)i + (long)j)*(long)this->dims_[2] + (long)k;
  }

  // Set the grid pointers for field
  /********************************************************************/
  void Field::assignGrid( double *x, double *y, double *z )
  /********************************************************************/
  {  
    this->x_ = x;
    this->y_ = y;
    this->z_ = z;
  }
  // Initialize the finite difference derivatives class
  /********************************************************************/
  void Field::derivFDInit( int order )
  /********************************************************************/
  {

    if( this->fd_ != NULL ) {
      cout << "finite difference class instance fd already allocated" << endl;
      exit (EXIT_FAILURE);
    }

    // First initialize the fd class instance
    this->fd_ = new FiniteDiff( this->dims_[0], this->x_, 
				this->dims_[1], this->y_, 
				this->dims_[2], this->z_, 
				order );
    
  }

  //////////////////////////////////////////////////////////////////////
  /// NON-CLASS MEMBER FUNCTIONS
  //////////////////////////////////////////////////////////////////////

}


