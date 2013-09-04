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
      cout << "Finite difference stencil set incorrectly\n" << endl;
      exit (EXIT_FAILURE);
    }

    for( int m=0; m<_fd.ssize; m++ ) {
      // Compute the distance from the 
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

    // First order
    if( _order == 1 ) {

      for (n=0; n<_ns; n++) {

	fdInit( 2, fd[n] );
	if( n== 0 ) {
	  fdSetStencil( 0, fd[n] ); // Forward differencing
	} else { 
	  fdSetStencil( -1, fd[n] ); // Backward differencing
	}
	fdSetDS( n, _s, fd[n] );
	fdSetCoef( 1, fd[n] );

      }

      // Second order 
    } else if( _order == 2 ) {
      
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
	
    } else if( _order == 3 ) {

	for (n=0; n<_ns; n++ ) {

	  fdInit( 4, fd[n] );

	  // Set the stencil
	  if( n == 0 ) {
	    // First point
	    fdSetStencil( 0, fd[n] ); // Forward differencing
	  } else if( n == 1 ) {
	    // Second point
	    fdSetStencil( -1, fd[n] ); // Forward-biased differencing
	  } else if ( n == _ns - 1 ) {
	    // Last point
	    fdSetStencil( -3, fd[n] ); // Backward differencing
	  } else {
	    // Interior points
	    fdSetStencil( -2, fd[n] ); // Backward-biased differencing
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

      cout << "Finite difference order >4 not supported" << endl;
      exit (EXIT_FAILURE);

    }
	

    return fd;
  }

}
