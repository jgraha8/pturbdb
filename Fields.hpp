#include "Derivs.hpp"

#ifndef FIELDS_H
#define FIELDS_H

namespace Fields {
  
  class Field {

  public:
    int ndim;
    int *dims;
    double *data;

    // FiniteDifference Class from the Derivs namespace
    Derivs::FiniteDiff *fd;
 
  public:

    // Constructor
    Field(){};
    Field( int _ndim, int *_dims );
    Field( const Field &g );
    // Deconstructor
    ~Field();

    // Assignment operator
    Field &operator=( const Field &a );
    Field &operator+=( const Field &a );
    Field &operator-=( const Field &a );
    Field &operator*=( const Field &a );
    Field &operator/=( const Field &a );

    // Get the size of the field
    long getSize();

    void add( const Field &a, const Field &b );
    void sub( const Field &a, const Field &b );
    void mul( const Field &a, const Field &b );
    void div( const Field &a, const Field &b );

    long index( int i, int j );
    long index( int i, int j, int k );

    // Derivative functions
    void derivFDInit( int _nx, double *_x, int _ny, double *_y, int _nz, double *_z, int _order );
    
  protected:

    void FieldInit( int *_dims );
    void FieldInit( int _ndim, int *_dims );

  }; 

}

#endif
