#include <iostream>
#include <cstdlib>
#include <string.h>
#include "Fields.hpp"

using namespace std;

namespace Fields {

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
  void Field::add( const Field &a, const Field &b )
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
  void Field::sub( const Field &a, const Field &b )
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
  void Field::mul( const Field &a, const Field &b ) 
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
  void Field::div( const Field &a, const Field &b ) 
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

  // Array index for 2D indexing
  /********************************************************************/
  long Field::index( int i, int j )
  /********************************************************************/
  {
    return (long)this->dims[1]*(long)i + (long)j;
  }
  // Array index for 3D indexing
  /********************************************************************/
  long Field::index( int i, int j, int k )
  /********************************************************************/
  {
    return ((long)this->dims[1]*(long)i + (long)j)*(long)this->dims[2] + (long)k;
  }

  // Initialize the finite difference derivatives class
  /********************************************************************/
  void Field::derivFDInit( int _nx, double *_x, int _ny, double *_y, int _nz, double *_z, int _order )
  /********************************************************************/
  {

    if( fd != NULL ) {
      cout << "finite difference class instance fd already allocated" << endl;
      exit (EXIT_FAILURE);
    }

    // First initialize the fd class instance
    fd = new Derivs::FiniteDiff();


    
  }

  //////////////////////////////////////////////////////////////////////
  /// NON-CLASS MEMBER FUNCTIONS
  //////////////////////////////////////////////////////////////////////

  // // Field Data Assignment
  // Field operator*( const Field &g, const Field &h )
  // {

  //   // Make a copy of class g
  //   Field f( g );

  //   cout << "f data addr " << f.data << endl;

  //   const long N=f.getSize();
  //   long indx=0;
  //   while( indx != N ) {
  //     f.data[indx] = g.data[indx] * h.data[indx];
  //     ++indx;
  //   }	
  //   return f;
  // }


}


