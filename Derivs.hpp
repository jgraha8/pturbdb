#include "autofd.h"

#ifndef DERIVS_H
#define DERIVS_H

// Extends the fields namespace
namespace Derivs 
{
  
  class Deriv 
  {
  public:
    // Constructor
    Deriv(){};
    // Deconstructor
    ~Deriv(){};
  };

  class FiniteDiff: public Deriv
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
    FiniteDiff( int _nx, double *_x, int _ny, double *_y, int _nz, double *_z, int _order );
    // Deconstructor
    ~FiniteDiff(){};

    void FiniteDiffInit( int _nx, double *_x, int _ny, double *_y, int _nz, double *_z, int _order );

  private:

    void fdInit( int _ssize, fd_t &_fd );
    void fdSetStencil( int _start_index, fd_t &_fd );
    void fdSetDS( int _n, double *_s, fd_t &_fd );
    void fdSetCoef( int _order_deriv, fd_t &_fd );
    fd_t *fdCreateD1( int _ns, double *_s, int _order );

  };


}

#endif
