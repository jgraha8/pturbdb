#include <iostream>
#include <cstdlib>
#include <string.h>

#include "mpi_field.hpp"

using namespace std;

namespace pturb_fields {
  
  //////////////////////////////////////////////////////////////////////
  /// CONSTRUCTORS
  //////////////////////////////////////////////////////////////////////
  /********************************************************************/
  MPIField::MPIField( MPI_Comm _comm, int _nproc, int _rank, 
		      MPIDecomp_t _decomp, int _ndim, const int *_dims )
  /********************************************************************/
  {
    MPIFieldInit(_comm, _nproc, _rank, _decomp, _ndim, _dims );
       
    // Get the processor local dimensions from the domain
    // decomposition; may have struct here for additional data
    int *_dims_loc = domainDecomp( _nproc, _rank, _decomp, _ndim, _dims );

    // Initialize the Field class
    FieldInit( _dims_loc );
    // Delete the local heap data
    delete [] _dims_loc;

  }

  /********************************************************************/
  MPIField::MPIField( const MPIField &g )
  /********************************************************************/
  {
    MPIFieldInit(g.comm, g.nproc, g.rank, g.decomp, g.ndim, g.dims_global);
    FieldInit( g.dims );
  }

  //////////////////////////////////////////////////////////////////////
  /// DECONSTRUCTOR
  //////////////////////////////////////////////////////////////////////

  /********************************************************************/
  MPIField::~MPIField(void) 
  /********************************************************************/
  {
    cout << "MPIField deconstructor" << endl;
    // Write the address of the field
    delete [] dims_global;
  }

  //////////////////////////////////////////////////////////////////////
  /// OPERATORS
  //////////////////////////////////////////////////////////////////////

  //////////////////////////////////////////////////////////////////////
  /// MEMBER FUNCTIONS
  //////////////////////////////////////////////////////////////////////

  /********************************************************************/
  long MPIField::getFieldAddr(void)
  /********************************************************************/
  {
    return (long)this;
  }

  /********************************************************************/
  long MPIField::getSizeGlobal(void)
  /********************************************************************/
  {
    long N=1;
    for (int n=0; n<ndim; n++ ) {
      N *= dims_global[n];
    }
    return N;
  }

  // void MPIField::setElementGlobal( double a, int i, int j, int k )
  // {
  //   int i_loc, j_loc, k_loc;
  //   indicesLoc( &i_loc, &j_loc, &k_loc );
  //   field->setElement( a, i_loc, j_loc, k_loc );
  // }

  /********************************************************************/
  void MPIField::MPIFieldInit(MPI_Comm _comm, int _nproc, int _rank, 
			      MPIDecomp_t _decomp, int _ndim, const int *_dims)
  /********************************************************************/
  {
    // Initialize data members
    comm   = _comm;
    nproc  = _nproc;
    rank   = _rank;
    decomp = _decomp;
    ndim   = _ndim; // Setting data in Field class

    dims_global = new int [ndim];   
    memcpy( dims_global, _dims, sizeof(*dims_global) * ndim );   
  }

  /********************************************************************/
  int *MPIField::domainDecomp( int _nproc, int _rank, MPIDecomp_t _decomp, 
			       int _ndim, const int *_dims )
  /********************************************************************/
  {

    int *dims_loc = new int[_ndim];

    // Now need to decompose 
    if( _decomp == SLAB_DECOMP ) {
      // Cut along the first dimension
      dims_loc[0] = _dims[0] / _nproc;
      for (int n=1; n<_ndim; n++ ) dims_loc[n] = _dims[n];
    } else {
      cout << "Only MPIDecomp_t SLAB_DECOMP currently supported" << endl;
      exit (EXIT_FAILURE);
    }    
    
    return dims_loc;

  }

  ////////////////////////////////////////////////////////////////////// 
  /// PRIVATE FUNCTIONS
  ////////////////////////////////////////////////////////////////////// 


  //////////////////////////////////////////////////////////////////////
  /// NON-CLASS MEMBER FUNCTIONS
  //////////////////////////////////////////////////////////////////////
  // // Field Data Assignment
  // MPIField operator*( const MPIField &a, const MPIField &b )
  // {

  //   // Make a copy of class a
  //   MPIField c( a );

  //   *c.field = *a.field * *b.field;
  //   cout << "c.field->data[1] " << c.field->data[1] << endl;
  //   return c;

  // }



  // MPIField MPIField::operator*( MPIField *g )
  // {
  //   // Create new MPIField
  //   //MPIField h( g->comm, g->nproc, g->rank, g->decomp, g->ndim, g->dims );
  //   return *this;
  //   //return h;
  // }

}


