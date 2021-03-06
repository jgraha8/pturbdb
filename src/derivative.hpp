#include "autofd.h"

#ifndef DERIVATIVE_H
#define DERIVATIVE_H

// Extends the fields namespace
namespace pturbdb
{
  
class Derivative
{

protected: 
	int nx_;
	int ny_;
	int nz_;

public:
	// Constructor
	Derivative(){};
	Derivative(int nx, int ny, int nz){ this->nx_ = nx; this->ny_ = ny; this->nz_ = nz; };
	// Deconstructor
	~Derivative(){};
};


class FiniteDiff: public Derivative
{

private:

	typedef struct {
		int n;
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
	FiniteDiff( int nx, const double *x, int ny, const double *y, int nz, const double *z, int order );
	// Deconstructor
	~FiniteDiff();

	void FiniteDiffInit( int nx, const double *x, int ny, const double *y, int nz, const double *z, int order );

	void ddx( int offset, int na, const double *a, int nda, double *da );
	void ddy( int offset, int na, const double *a, int nda, double *da );
	void ddz( int offset, int na, const double *a, int nda, double *da );
	void d2dx2( int offset, int na, const double *a, int nda, double *da );
	void d2dy2( int offset, int na, const double *a, int nda, double *da );
	void d2dz2( int offset, int na, const double *a, int nda, double *da );

private:

	void fdInit( int ssize, fd_t &fd );

	void fdSetStencil( int start_index, fd_t &fd );
	void fdSetDS( int n, const double *s, fd_t &fd );
	void fdSetCoef( int order_deriv, fd_t &fd );
	fd_t *fdCreateD1( int ns, const double *s, int order );
	fd_t *fdCreateD2( int ns, const double *s, int order );
	void fdDelete( int n, fd_t **fd );
	// Performs the finite differencing operation
	void fdOp( const fd_t *fd, int offset, int na, const double *a, int nda, double *da );

};


}

#endif
