#include <iostream>
#include <cstdlib>
#include <string.h>
#include "field.hpp"

using namespace std;

namespace pturb_fields {

    // Field class member functions
  /********************************************************************/
  Field::Field( int _ndim, int *_dims ) 
  /********************************************************************/
  {
    cout << "Initializing Field" << endl;
    FieldInit( _ndim, _dims );
  }

  /********************************************************************/
  Field::Field( const Field &g ) 
  /********************************************************************/
  {
    cout << "Copy Constructor Initializing Field" << endl;

    FieldInit( g.ndim, g.dims );
    // Copy the field data from g to this->data
    memcpy(this->data,g.data,sizeof(*this->data)*this->getSize());

  }

  /********************************************************************/
  Field::~Field(void) 
  /********************************************************************/
  {
    cout << "Field deconstructor" << endl;
    delete[] dims;
    delete[] data;
  }

  // Field class initializer; expects ndim to already be set.
  // Initilizes:
  //   dims
  //   data
  /********************************************************************/
  void Field::FieldInit( int *_dims ) 
  /********************************************************************/
  {
    dims = new int [ndim];    
    // Copy the dims vector
    memcpy( dims, _dims, sizeof(*dims) * ndim );

    // Initialize the data field
    long N = getSize();
    data = new double[N];

    // Initialize pointers to null
    fd = NULL;
  }

  // Field class initializer
  // Sets:
  //   ndim
  // Initilizes:
  //   dims
  //   data
  /********************************************************************/
  void Field::FieldInit( int _ndim, int *_dims ) 
  /********************************************************************/
  {
    ndim = _ndim;
    FieldInit( _dims );
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
      memcpy( data, a.data, sizeof(*data) * getSize() );
    }	
    return *this;
  }


  /********************************************************************/
  Field& Field::operator+=( const Field &a )
  /********************************************************************/
  {
    const long N=getSize();
    long indx=0;
    while( indx != N ) {
      data[indx] += a.data[indx];
      ++indx;
    }	
    return *this;
  }

  /********************************************************************/
  Field& Field::operator-=( const Field &a )
  /********************************************************************/
  {
    const long N=getSize();
    long indx=0;
    while( indx != N ) {
      data[indx] -= a.data[indx];
      ++indx;
    }	
    return *this;
  }

  /********************************************************************/
  Field& Field::operator*=( const Field &a )
  /********************************************************************/
  {
    const long N=getSize();
    long indx=0;
    while( indx != N ) {
      data[indx] *= a.data[indx];
      ++indx;
    }	
    return *this;
  }

  /********************************************************************/
  Field& Field::operator/=( const Field &a )
  /********************************************************************/
  {
    const long N=getSize();
    long indx=0;
    while( indx != N ) {
      data[indx] /= a.data[indx];
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
    for (int n=0; n<ndim; n++ ) {
      N *= dims[n];
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
      this->data[indx] = a.data[indx] + b.data[indx];
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
      this->data[indx] = a.data[indx] - b.data[indx];
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
      this->data[indx] = a.data[indx] * b.data[indx];
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
      this->data[indx] = a.data[indx] / b.data[indx];
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
    dndxn( dd, a );
  }

  /********************************************************************/
  void Field::ddy( Field &a ) 
  /********************************************************************/
  {
    void (FiniteDiff::*dd)(int, double *, double *);
    dd = &FiniteDiff::ddy;
    dndyn( dd, a );
  }

  /********************************************************************/
  void Field::ddz( Field &a ) 
  /********************************************************************/
  {
    void (FiniteDiff::*dd)(int, double *, double *);
    dd = &FiniteDiff::ddz;
    dndzn( dd, a );
  }

  /********************************************************************/
  void Field::d2dx2( Field &a ) 
  /********************************************************************/
  {
    void (FiniteDiff::*dd)(int, double *, double *);
    dd = &FiniteDiff::d2dx2;
    dndxn( dd, a );
  }
  /********************************************************************/
  void Field::d2dy2( Field &a ) 
  /********************************************************************/
  {
    void (FiniteDiff::*dd)(int, double *, double *);
    dd = &FiniteDiff::d2dy2;
    dndyn( dd, a );
  }

  /********************************************************************/
  void Field::d2dz2( Field &a ) 
  /********************************************************************/
  {
    void (FiniteDiff::*dd)(int, double *, double *);
    dd = &FiniteDiff::d2dz2;
    dndzn( dd, a );
  }

  /********************************************************************/
  void Field::d2dxy( Field &a ) 
  /********************************************************************/
  {
    void (FiniteDiff::*dd)(int, double *, double *);
    dd = &FiniteDiff::ddx;
    dndxn( dd, a );
    dd = &FiniteDiff::ddy;
    dndyn( dd, *this );
  } 
  /********************************************************************/
  void Field::d2dxz( Field &a ) 
  /********************************************************************/
  {
    void (FiniteDiff::*dd)(int, double *, double *);
    dd = &FiniteDiff::ddx;
    dndxn( dd, a );
    dd = &FiniteDiff::ddz;
    dndzn( dd, *this );
  }
  /********************************************************************/
  void Field::d2dyz( Field &a ) 
  /********************************************************************/
  {
    void (FiniteDiff::*dd)(int, double *, double *);
    dd = &FiniteDiff::ddy;
    dndyn( dd, a );
    dd = &FiniteDiff::ddz;
    dndzn( dd, *this );
  }    

  /********************************************************************/
  void Field::dndxn( void (FiniteDiff::*dd)( int, double *, double *), Field &a ) 
  /********************************************************************/
  {

    // Create buffer to store x-data
    const int nx = this->dims[0];
    const int ny = this->dims[1];
    const int nz = this->dims[2];
    
#ifdef BOUNDS_CHECK
    // First check that the fields are the same size
    if ( nx != a.dims[0] || ny != a.dims[1] || nz != a.dims[2] ) {
      cout << "Mismatch in field sizes" << endl;
      exit (EXIT_FAILURE);
    }   
#endif

    double *ax = new double[nx];
    double *dax = new double[nx];
    
    for (int j=0; j<ny; j++ ) {
      for (int k=0; k<nz; k++ ) {
	// Pack buffer
	for ( int i=0; i<nx; i++ ) ax[i] = a.data[a.index(i,j,k)];
	// Compute the derivative using the FiniteDiff class
	(this->fd->*dd)( nx, ax, dax );
	// Unpack buffer
	for ( int i=0; i<nx; i++ ) this->data[this->index(i,j,k)] = dax[i];
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
    const int nx = this->dims[0];
    const int ny = this->dims[1];
    const int nz = this->dims[2];
    
#ifdef BOUNDS_CHECK
    // First check that the fields are the same size
    if ( nx != a.dims[0] || ny != a.dims[1] || nz != a.dims[2] ) {
      cout << "Mismatch in field sizes" << endl;
      exit (EXIT_FAILURE);
    }   
#endif

    double *ay = new double[ny];
    double *day = new double[ny];
    
    for ( int i=0; i<nx; i++ ) {
      for (int k=0; k<nz; k++ ) {
	// Pack buffer
	for (int j=0; j<ny; j++ ) ay[j] = a.data[a.index(i,j,k)];
	// Compute the derivative using the FiniteDiff class
	(this->fd->*dd)( ny, ay, day );
	// Unpack buffer
	for ( int j=0; j<ny; j++ ) this->data[this->index(i,j,k)] = day[j];
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
    const int nx = this->dims[0];
    const int ny = this->dims[1];
    const int nz = this->dims[2];
    
#ifdef BOUNDS_CHECK
    // First check that the fields are the same size
    if ( nx != a.dims[0] || ny != a.dims[1] || nz != a.dims[2] ) {
      cout << "Mismatch in field sizes" << endl;
      exit (EXIT_FAILURE);
    }   
#endif

    double *az = new double[nz];
    double *daz = new double[nz];
    
    for ( int i=0; i<nx; i++ ) {
      for (int j=0; j<ny; j++ ) {
	// Pack buffer
	for (int k=0; k<nz; k++ ) az[k] = a.data[a.index(i,j,k)];
	// Compute the derivative using the FiniteDiff class
	(this->fd->*dd)( nz, az, daz );
	// Unpack buffer
	for (int k=0; k<nz; k++ ) this->data[this->index(i,j,k)] = daz[k];
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
    return ((long)dims[1]*(long)i + (long)j)*(long)dims[2] + (long)k;
  }

  // Set the grid pointers for field
  /********************************************************************/
  void Field::assignGrid( double *_x, double *_y, double *_z )
  /********************************************************************/
  {  
    x = _x;
    y = _y;
    z = _z;
  }
  // Initialize the finite difference derivatives class
  /********************************************************************/
  void Field::derivFDInit( int _order )
  /********************************************************************/
  {

    if( fd != NULL ) {
      cout << "finite difference class instance fd already allocated" << endl;
      exit (EXIT_FAILURE);
    }

    // First initialize the fd class instance
    fd = new FiniteDiff( dims[0], x, dims[1], y, dims[2], z, _order );
    
  }

  //////////////////////////////////////////////////////////////////////
  /// NON-CLASS MEMBER FUNCTIONS
  //////////////////////////////////////////////////////////////////////

}


