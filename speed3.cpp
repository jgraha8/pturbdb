#include <iostream>
//#include <stdlib.h>
#include <stdlib.h>

#define BOUNDS_CHECK

#include "MPIFields.hpp"

#define NX 512
#define NY 512
#define NZ 512

using namespace std;
using namespace MPIFields;

int main ( int argc, char *argv[]) 
{
  int id;
  int ierr;
  int p;
  double wtime;

  ierr = MPI_Init ( &argc, &argv );
  ierr = MPI_Comm_size ( MPI_COMM_WORLD, &p );
  ierr = MPI_Comm_rank ( MPI_COMM_WORLD, &id );

  int ndim=3;
  int dims[]={ NX, NY, NZ }; 

  MPIField *f = new MPIField( MPI_COMM_WORLD, p, id, SLAB_DECOMP, ndim, dims );
  MPIField *g = new MPIField( *f ); // Make a copy of f
				   
  for (int n=0; n<p; n++) {
    if( n == id ) {
      cout << "Number of local data elements : " << g->getSize() << endl;
      cout << "Number of global data elements : " << g->getSizeGlobal() << endl;
      cout << "f field addr : " << g->getFieldAddr() << endl;
    }
    ierr = MPI_Barrier( MPI_COMM_WORLD );
  }

  for (int i=0; i<f->dims[0]; i++) {
    for (int j=0; j<f->dims[1]; j++) {
      for (int k=0; k<f->dims[2]; k++) {
	f->data[f->index(i,j,k)] = (double)i*j + k;
      }
    }
  }
  
  double *x = new double[f->dims[0]];
  double *y = new double[f->dims[1]];
  double *z = new double[f->dims[2]];

  for (int i=0; i<f->dims[0]; i++ ) x[i] = (double)i;
  for (int j=0; j<f->dims[0]; j++ ) y[j] = (double)j;
  for (int k=0; k<f->dims[0]; k++ ) z[k] = (double)k;

  // Initialize derivatives for f
  f->derivFDInit( f->dims[0], x, f->dims[1], y, f->dims[2], z, 4 );

  // Print the FD coefficients along the x, y, and z directions
  for (int i=0; i<f->dims[0]; i++ ) {
    cout << "i : " << i;
    for (int m=0; m<f->fd->fd_ddx[i].ssize; m++ ) {
      cout << " " << f->fd->fd_ddx[i].coef[m];
    }
    cout << endl;
  }


  // Set g to have same field data (does not copy entire class;
  // assumes they are initilized the same);
  (*g)=(*f);

  // Now take product of f and g and assign to f
  f->add(*f,*g);
  f->sub(*f,*g);
  f->mul(*f,*g);
  f->div(*f,*g); 

  if ( id == 0 ) 
    {
      wtime = MPI_Wtime ( );

      printf ( "\n" );
      printf ( "HELLO_MPI - Master process:\n" );
      printf ( "  C/MPI version\n" );
      printf ( "  An MPI example program.\n" );
      printf ( "\n" );
      printf ( "  The number of processes is %d.\n", p );
      printf ( "\n" );
    }
  /*
  Every process prints a hello.
  */
  printf ( "  Process %d says 'Hello, world!'\n", id );

  if ( id == 0 )
    {
      wtime = MPI_Wtime ( ) - wtime;
      printf ( "  Elapsed wall clock time = %f seconds.\n", wtime );
    }
  
  // Terminate MPI.
   ierr = MPI_Finalize ( );
  
  if ( id == 0 )
    {
      printf ( "\n" );
      printf ( "HELLO_MPI - Master process:\n" );
      printf ( "  Normal end of execution: 'Goodbye, world!'\n" );
      printf ( "\n" );
    }

  return 0;
}
