#include <iostream>
#include <cstdlib>
#include <string.h>
#include "Fields.hpp"

using namespace std;

namespace Derivs {

  FiniteDiff::FiniteDiff( int _nx, double *_x, int _ny, double *_y, int _nz, double *_z, int _order ) 
  {
    FiniteDiffInit( _nx, _x, _ny, _y, _nz, _z, _order );
  }    

  void FiniteDiff::FiniteDiffInit( int _nx, double *_x, int _ny, double *_y, int _nz, double *_z, int _order )
  {

    fd_ddx = fdCreateD1( _nx, _x, _order );
    fd_ddy = fdCreateD1( _ny, _y, _order );
    fd_ddz = fdCreateD1( _nz, _z, _order );
    
    fd_d2dx2 = fdCreateD2( _nx, _x, _order );
    fd_d2dy2 = fdCreateD2( _ny, _y, _order );
    fd_d2dz2 = fdCreateD2( _nz, _z, _order );

  }

  /********************************************************************/
  void FiniteDiff::fdInit( int _ssize, fd_t &_fd )
  /********************************************************************/
  {
    // Allocate the members
    _fd.ssize   = _ssize;
    _fd.stencil = new int[_fd.ssize];
    _fd.ds      = new autofd_real_t[_fd.ssize];
    _fd.coef    = new autofd_real_t[_fd.ssize];
  }

  /********************************************************************/
  void FiniteDiff::fdSetStencil( int _start_index, fd_t &_fd )
  /********************************************************************/
  {
    for ( int m=0; m<_fd.ssize; m++ ) _fd.stencil[m] = _start_index + m;
  }

  /********************************************************************/
  void FiniteDiff::fdSetDS( int _n, double *_s, fd_t &_fd )
  /********************************************************************/
  {

    int zero_index=-1;
    // First find the zero index of the stencil
    for( int m=0; m<_fd.ssize; m++ ) {
      if( _fd.stencil[m] == 0 ) {
	zero_index = m;
	break;
      }
    }

    if( zero_index == -1 ) {
      // Did not find the zero index 
      cout << "Finite difference stencil set incorrectly\n" << endl;
      exit (EXIT_FAILURE);
    }

    for( int m=0; m<_fd.ssize; m++ ) {
      // Compute the distance from the zero index location
      _fd.ds[m] = _s[ _n + _fd.stencil[m] ] - _s[ _n + _fd.stencil[zero_index] ];
    }

  }

  /********************************************************************/
  void FiniteDiff::fdSetCoef( int _order_deriv, fd_t &_fd )
  /********************************************************************/
  {
    // Generate the finite difference coefficients
    autofd_stencil_coefs_c( _fd.ssize, _fd.stencil, _fd.ds, _order_deriv, _fd.coef );
  }

  /********************************************************************/
  FiniteDiff::fd_t *FiniteDiff::fdCreateD1( int _ns, double *_s, int _order ) 
  /********************************************************************/
  {

    FiniteDiff::fd_t *fd = new FiniteDiff::fd_t[_ns];

    int n=0;

    // Second order 
    if( _order == 2 ) {
      
	for (n=0; n<_ns; n++ ) {

	  fdInit( 3, fd[n] );

	  // Set the stencil
	  if( n == 0 ) {
	    // First point
	    fdSetStencil( 0, fd[n] ); // Forward differencing
	  } else if ( n == _ns - 1 ) {
	    // Last point
	    fdSetStencil( -2, fd[n] ); // Backward differencing
	  } else {
	    // Interior points
	    fdSetStencil( -1, fd[n] ); // Central differencing
	  }
	  fdSetDS( n, _s, fd[n] );
	  fdSetCoef( 1, fd[n] );

	}
	
    } else if( _order == 4 ) {

	for (n=0; n<_ns; n++ ) {

	  fdInit( 5, fd[n] );

	  // Set the stencil
	  if( n == 0 ) {
	    // First point
	    fdSetStencil( 0, fd[n] ); // Forward differencing
	  } else if( n == 1 ) {
	    // Second point
	    fdSetStencil( -1, fd[n] ); // Forward-biased differencing
	  } else if ( n == _ns - 2 ) {
	    // Second to last point
	    fdSetStencil( -3, fd[n] ); // Backward-biased differencing
	  } else if ( n == _ns - 1 ) {
	    // Last point
	    fdSetStencil( -4, fd[n] ); // Backward differencing
	  } else {
	    // Interior points
	    fdSetStencil( -2, fd[n] ); // Central differencing
	  }

	  fdSetDS( n, _s, fd[n] );
	  fdSetCoef( 1, fd[n] );
	  
	}

    } else {

      cout << "Finite difference order 2 and 4 only supported" << endl;
      exit (EXIT_FAILURE);

    }
    return fd;
  }

  /********************************************************************/
  FiniteDiff::fd_t *FiniteDiff::fdCreateD2( int _ns, double *_s, int _order ) 
  /********************************************************************/
  {

    FiniteDiff::fd_t *fd = new FiniteDiff::fd_t[_ns];

    int n=0;

    // Second order 
    if( _order == 2 ) {
      
	for (n=0; n<_ns; n++ ) {

	  // Set the stencil
	  if( n == 0 ) {
	    // First point
	    fdInit( 4, fd[n] );
	    fdSetStencil( 0, fd[n] ); // Forward differencing
	  } else if ( n == _ns - 1 ) {
	    // Last point
	    fdInit( 4, fd[n] );
	    fdSetStencil( -3, fd[n] ); // Backward differencing
	  } else {
	    // Interior points
	    fdInit( 3, fd[n] );
	    fdSetStencil( -1, fd[n] ); // Central differencing
	  }
	  fdSetDS( n, _s, fd[n] );
	  fdSetCoef( 2, fd[n] );

	}
	
    } else if( _order == 4 ) {

	for (n=0; n<_ns; n++ ) {

	  // Set the stencil
	  if( n == 0 ) {
	    // First point
	    fdInit( 6, fd[n] );
	    fdSetStencil( 0, fd[n] ); // Forward differencing
	  } else if( n == 1 ) {
	    // Second point
	    fdInit( 6, fd[n] );
	    fdSetStencil( -1, fd[n] ); // Forward-biased differencing
	  } else if ( n == _ns - 2 ) {
	    // Second to last point
	    fdInit( 6, fd[n] );
	    fdSetStencil( -4, fd[n] ); // Backward-biased differencing
	  } else if ( n == _ns - 1 ) {
	    // Last point
	    fdInit( 6, fd[n] );
	    fdSetStencil( -5, fd[n] ); // Backward differencing
	  } else {
	    // Interior points
	    fdInit( 5, fd[n] );
	    fdSetStencil( -2, fd[n] ); // Central differencing
	  }

	  fdSetDS( n, _s, fd[n] );
	  fdSetCoef( 2, fd[n] );
	  
	}

    } else {

      cout << "Finite difference order 2 and 4 only are supported" << endl;
      exit (EXIT_FAILURE);

    }
	
    return fd;
  }

  //////////////////////////////////////////////////////////////////////
  /// DERIVATIVES
  ////////////////////////////////////////////////////////////////////// 

  /********************************************************************/
  void FiniteDiff::ddx( int _n, double *_a, double *_da )
  /********************************************************************/
  {
    this->fdOp( this->fd_ddx, _n, _a, _da );
  }
  /********************************************************************/
  void FiniteDiff::ddy( int _n, double *_a, double *_da )
  /********************************************************************/
  {
    this->fdOp( this->fd_ddy, _n, _a, _da );
  }
  /********************************************************************/
  void FiniteDiff::ddz( int _n, double *_a, double *_da )
  /********************************************************************/
  {
    this->fdOp( this->fd_ddz, _n, _a, _da );
  }
  /********************************************************************/
  void FiniteDiff::d2dx2( int _n, double *_a, double *_da )
  /********************************************************************/
  {
    this->fdOp( this->fd_d2dx2, _n, _a, _da );
  }
  /********************************************************************/
  void FiniteDiff::d2dy2( int _n, double *_a, double *_da )
  /********************************************************************/
  {
    this->fdOp( this->fd_d2dy2, _n, _a, _da );
  }
  /********************************************************************/
  void FiniteDiff::d2dz2( int _n, double *_a, double *_da )
  /********************************************************************/
  {
    this->fdOp( this->fd_d2dz2, _n, _a, _da );
  }

  //
  // Finite difference operation function
  //
  /********************************************************************/
  void FiniteDiff::fdOp( FiniteDiff::fd_t *fd, int _n, double *_a, double *_da )
  /********************************************************************/
  {
    for ( int i=0; i<_n; i++ ) {
      _da[i] = autofd_derivative_c( fd[i].ssize, fd[i].stencil,  fd[i].coef,  _n, _a, i);
    }
  } 

}
