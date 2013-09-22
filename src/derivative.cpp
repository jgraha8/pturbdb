#include <iostream>
#include <cstdlib>
#include <cstdio>
#include <string.h>
#include "derivative.hpp"

using namespace std;

namespace pturbdb {

FiniteDiff::FiniteDiff( int nx, double *x, int ny, double *y, int nz, double *z, int order ) 
{
	this->fd_ddx = NULL;
	this->fd_ddy = NULL;
	this->fd_ddz = NULL;
	
	this->fd_d2dx2 = NULL;
	this->fd_d2dy2 = NULL;
	this->fd_d2dz2 = NULL;
	
	FiniteDiffInit( nx, x, ny, y, nz, z, order );
}    

void FiniteDiff::FiniteDiffInit( int nx, double *x, int ny, double *y, int nz, double *z, int order )
{
#ifdef VERBOSE
	printf("FiniteDiff::FiniteDiffInit: entering\n");
#endif

	this->fd_ddx = this->fdCreateD1( nx, x, order );
	this->fd_ddy = this->fdCreateD1( ny, y, order );
	this->fd_ddz = this->fdCreateD1( nz, z, order );
	  
	this->fd_d2dx2 = this->fdCreateD2( nx, x, order );
	this->fd_d2dy2 = this->fdCreateD2( ny, y, order );
	this->fd_d2dz2 = this->fdCreateD2( nz, z, order );

#ifdef VERBOSE
	printf("FiniteDiff::FiniteDiffInit: exiting\n");
#endif

}

/********************************************************************/
void FiniteDiff::fdInit( int ssize, fd_t &fd )
/********************************************************************/
{
	// Allocate the members
	fd.ssize   = ssize;
	fd.stencil = new int[fd.ssize];
	fd.ds      = new autofd_real_t[fd.ssize];
	fd.coef    = new autofd_real_t[fd.ssize];
}

/********************************************************************/
void FiniteDiff::fdSetStencil( int start_index, fd_t &fd )
/********************************************************************/
{
	for ( int m=0; m<fd.ssize; m++ ) fd.stencil[m] = start_index + m;
}

/********************************************************************/
void FiniteDiff::fdSetDS( int n, double *s, fd_t &fd )
/********************************************************************/
{

	int zero_index=-1;
	// First find the zero index of the stencil
	for( int m=0; m<fd.ssize; m++ ) {
		if( fd.stencil[m] == 0 ) {
			zero_index = m;
			break;
		}
	}

	if( zero_index == -1 ) {
		// Did not find the zero index 
		cout << "Finite difference stencil set incorrectly\n" << endl;
		exit (EXIT_FAILURE);
	}

	for( int m=0; m<fd.ssize; m++ ) {
		// Compute the distance from the zero index location
		fd.ds[m] = s[ n + fd.stencil[m] ] - s[ n + fd.stencil[zero_index] ];
	}

}

/********************************************************************/
void FiniteDiff::fdSetCoef( int order_deriv, fd_t &fd )
/********************************************************************/
{
	// Generate the finite difference coefficients
	autofd_stencil_coefs_c( fd.ssize, fd.stencil, fd.ds, order_deriv, fd.coef );
}

/********************************************************************/
FiniteDiff::fd_t *FiniteDiff::fdCreateD1( int ns, double *s, int order ) 
/********************************************************************/
{

	FiniteDiff::fd_t *fd = new FiniteDiff::fd_t[ns];

	int n=0;

	// Second order 
	if( order == 2 ) {
      
		for (n=0; n<ns; n++ ) {

			fdInit( 3, fd[n] );

			// Set the stencil
			if( n == 0 ) {
				// First point
				fdSetStencil( 0, fd[n] ); // Forward differencing
			} else if ( n == ns - 1 ) {
				// Last point
				fdSetStencil( -2, fd[n] ); // Backward differencing
			} else {
				// Interior points
				fdSetStencil( -1, fd[n] ); // Central differencing
			}
			fdSetDS( n, s, fd[n] );
			fdSetCoef( 1, fd[n] );

		}
	
	} else if( order == 4 ) {

		for (n=0; n<ns; n++ ) {

			fdInit( 5, fd[n] );

			// Set the stencil
			if( n == 0 ) {
				// First point
				fdSetStencil( 0, fd[n] ); // Forward differencing
			} else if( n == 1 ) {
				// Second point
				fdSetStencil( -1, fd[n] ); // Forward-biased differencing
			} else if ( n == ns - 2 ) {
				// Second to last point
				fdSetStencil( -3, fd[n] ); // Backward-biased differencing
			} else if ( n == ns - 1 ) {
				// Last point
				fdSetStencil( -4, fd[n] ); // Backward differencing
			} else {
				// Interior points
				fdSetStencil( -2, fd[n] ); // Central differencing
			}

			fdSetDS( n, s, fd[n] );
			fdSetCoef( 1, fd[n] );
	  
		}

	} else {

		cout << "Finite difference order 2 and 4 only supported" << endl;
		exit (EXIT_FAILURE);

	}
	return fd;
}

/********************************************************************/
FiniteDiff::fd_t *FiniteDiff::fdCreateD2( int ns, double *s, int order ) 
/********************************************************************/
{

	FiniteDiff::fd_t *fd = new FiniteDiff::fd_t[ns];

	int n=0;

	// Second order 
	if( order == 2 ) {
      
		for (n=0; n<ns; n++ ) {

			// Set the stencil
			if( n == 0 ) {
				// First point
				fdInit( 4, fd[n] );
				fdSetStencil( 0, fd[n] ); // Forward differencing
			} else if ( n == ns - 1 ) {
				// Last point
				fdInit( 4, fd[n] );
				fdSetStencil( -3, fd[n] ); // Backward differencing
			} else {
				// Interior points
				fdInit( 3, fd[n] );
				fdSetStencil( -1, fd[n] ); // Central differencing
			}
			fdSetDS( n, s, fd[n] );
			fdSetCoef( 2, fd[n] );

		}
	
	} else if( order == 4 ) {

		for (n=0; n<ns; n++ ) {

			// Set the stencil
			if( n == 0 ) {
				// First point
				fdInit( 6, fd[n] );
				fdSetStencil( 0, fd[n] ); // Forward differencing
			} else if( n == 1 ) {
				// Second point
				fdInit( 6, fd[n] );
				fdSetStencil( -1, fd[n] ); // Forward-biased differencing
			} else if ( n == ns - 2 ) {
				// Second to last point
				fdInit( 6, fd[n] );
				fdSetStencil( -4, fd[n] ); // Backward-biased differencing
			} else if ( n == ns - 1 ) {
				// Last point
				fdInit( 6, fd[n] );
				fdSetStencil( -5, fd[n] ); // Backward differencing
			} else {
				// Interior points
				fdInit( 5, fd[n] );
				fdSetStencil( -2, fd[n] ); // Central differencing
			}

			fdSetDS( n, s, fd[n] );
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

// 
// The offset variable allows the data in a to be a sub-vector of the vector
// spanned by the finite differencing struct.
// 

/********************************************************************/
void FiniteDiff::ddx( int offset, int na, double *a, int nda, double *da )
/********************************************************************/
{
#ifdef VERBOSE
	printf("FiniteDiff::ddx: entering\n");
#endif
	if( this->fd_ddx == NULL ) {
		printf("FiniteDiff::fdOp: received null finite difference struct\n");
	}

	this->fdOp( this->fd_ddx, offset, na, a, nda, da );
#ifdef VERBOSE
	printf("FiniteDiff::ddx: exiting\n");
#endif
}
/********************************************************************/
void FiniteDiff::ddy( int offset, int na, double *a, int nda, double *da )
/********************************************************************/
{
	this->fdOp( this->fd_ddy, offset, na, a, nda, da );
}
/********************************************************************/
void FiniteDiff::ddz( int offset, int na, double *a, int nda, double *da )
/********************************************************************/
{
	this->fdOp( this->fd_ddz, offset, na, a, nda, da );
}
/********************************************************************/
void FiniteDiff::d2dx2( int offset, int na, double *a, int nda, double *da )
/********************************************************************/
{
	this->fdOp( this->fd_d2dx2, offset, na, a, nda, da );
}
/********************************************************************/
void FiniteDiff::d2dy2( int offset, int na, double *a, int nda, double *da )
/********************************************************************/
{
	this->fdOp( this->fd_d2dy2, offset, na, a, nda, da );
}
/********************************************************************/
void FiniteDiff::d2dz2( int offset, int na, double *a, int nda, double *da )
/********************************************************************/
{
	this->fdOp( this->fd_d2dz2, offset, na, a, nda, da );
}

//
// Finite difference operation function
//
/********************************************************************/
void FiniteDiff::fdOp( FiniteDiff::fd_t *fd, int offset, int na, double *a, int nda, double *da )
/********************************************************************/
{
#ifdef VERBOSE
	printf("FiniteDiff::fdOp: entering\n");
#ifdef DEBUG
	printf("    offset = %d, na = %d, nda = %d\n", offset, na, nda);
#endif
#endif
	if( fd == NULL ) {
		printf("FiniteDiff::fdOp: received null finite difference struct\n");
	}
	for ( int i=0; i<nda; i++ ) {
		da[i] = autofd_derivative_c( fd[i+offset].ssize, 
					     fd[i+offset].stencil,  
					     fd[i+offset].coef,  
					     na, a, 
					     i+offset);
	}
#ifdef VERBOSE
	printf("FiniteDiff::fdOp: exiting\n");
#endif


} 

}
