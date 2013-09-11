#include "autofd.h"

#ifndef DERIVATIVE_H
#define DERIVATIVE_H

// Extends the fields namespace
namespace pturb_fields
{
  
  class Derivative
  {

  public:
    // Constructor
    Derivative(){};
    // Deconstructor
    ~Derivative(){};
  };


  class FiniteDiff: public Derivative
  {

  private:

    typedef struct {
      int ssize;
      int *stencil; 
      autofd_real_t *ds;
      autofd_real_t *coef;  // Coefficients
    } fd_t;

    fd_t *fd_ddx;
    fd_t *fd_ddy;
    fd_t *fd_ddz;

    fd_t *fd_d2dx2;
    fd_t *fd_d2dy2;
    fd_t *fd_d2dz2;
    
  public:

    // Constructor
    FiniteDiff(){};
    FiniteDiff( int nx, double *x, int ny, double *y, int nz, double *z, int order );
    // Deconstructor
    ~FiniteDiff(){};

    void FiniteDiffInit( int nx, double *x, int ny, double *y, int nz, double *z, int order );

    void ddx( int offset, int n, double *a, double *da );
    void ddy( int offset, int n, double *a, double *da );
    void ddz( int offset, int n, double *a, double *da );
    void d2dx2( int offset, int n, double *a, double *da );
    void d2dy2( int offset, int n, double *a, double *da );
    void d2dz2( int offset, int n, double *a, double *da );

  private:

    void fdInit( int ssize, fd_t &fd );
    void fdSetStencil( int start_index, fd_t &fd );
    void fdSetDS( int n, double *s, fd_t &fd );
    void fdSetCoef( int order_deriv, fd_t &fd );
    fd_t *fdCreateD1( int ns, double *s, int order );
    fd_t *fdCreateD2( int ns, double *s, int order );

    // Performs the finite differencing operation
    void fdOp( fd_t *fd, int offset, int n, double *a, double *da );

  };


}

#endif
