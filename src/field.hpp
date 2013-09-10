#include "derivative.hpp"

#ifndef FIELD_H
#define FIELD_H

namespace pturb_fields {
  
  class Field {

  private:
    int ndim_;
    int *dims_;
    double *data_;

    // Grid pointers; these are not allocated but must be set with setGrid
    double *x_, *y_, *z_;

    // FiniteDiff Class from the D namespace
    FiniteDiff *fd_;
 
  public:

    // Constructor
    Field(){};
    Field( int ndim, int *dims );
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

    void add( Field &a, Field &b );
    void sub( Field &a, Field &b );
    void mul( Field &a, Field &b );
    void div( Field &a, Field &b );

    void ddx( Field &a );
    void ddy( Field &a );
    void ddz( Field &a );

    void d2dx2( Field &a );
    void d2dy2( Field &a );
    void d2dz2( Field &a );
    void d2dxy( Field &a );
    void d2dxz( Field &a );
    void d2dyz( Field &a );

    long index( int i, int j, int k );

    // Grid 
    void assignGrid( double *x, double *y, double *z );
    // Derivative functions
    void derivFDInit( int order );
    
  protected:

    void FieldInit( int *dims );
    void FieldInit( int ndim, int *dims );

    void dndxn( void (FiniteDiff::*dd)( int, double *, double *), Field &a );
    void dndyn( void (FiniteDiff::*dd)( int, double *, double *), Field &a );
    void dndzn( void (FiniteDiff::*dd)( int, double *, double *), Field &a );

  }; 

}

#endif
