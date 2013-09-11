#include <iostream>
//#include <stdlib.h>
#include <stdlib.h>

#include "field.hpp"

#define NX 128
#define NY 128
#define NZ 128

using namespace std;
using namespace pturb_fields;

int main ( int argc, char *argv[]) 
{
  double wtime;

  MPI_Init ( &argc, &argv );
  //  ierr = MPI_Comm_size ( MPI_COMM_WORLD, &p );
  //ierr = MPI_Comm_rank ( MPI_COMM_WORLD, &id );

  int ndim=3;
  int dims[]={ NX, NY, NZ }; 
  int periodic[] = {1, 1, 1};

  Field *f = new Field( dims, FIELD_DECOMP_PENCIL, periodic, 2 );
  Field *g = new Field( *f ); // Make a copy of f

  int rank = f->getMpiTopology()->rank;
  int nproc = f->getMpiTopology()->nproc;
				   
  for (int n=0; n<nproc; n++) {
    if( n == rank ) {
      cout << "Number of global data elements : " << g->getSize() << endl;
      cout << "Number of local data elements : " << g->getSizeLocal() << endl;
      cout << "Number of operation data elements : " << g->getSizeOperation() << endl;
    }
    MPI_Barrier( g->getMpiTopology()->comm );
  }

 // for (int i=0; i<f->dims[0]; i++) {
 //    for (int j=0; j<f->dims[1]; j++) {
 //      for (int k=0; k<f->dims[2]; k++) {
 // 	f->data[f->index(i,j,k)] = (double)i*j + k;
 //      }
 //    }
 //  }
  
 //  double *x = new double[f->dims[0]];
 //  double *y = new double[f->dims[1]];
 //  double *z = new double[f->dims[2]];

 //  for (int i=0; i<f->dims[0]; i++ ) x[i] = (double)i;
 //  for (int j=0; j<f->dims[0]; j++ ) y[j] = (double)j;
 //  for (int k=0; k<f->dims[0]; k++ ) z[k] = (double)k;

 //  // Assign the grid pointers; needed for finite differencing
 //  f->assignGrid( x, y, z );

 //  // Initialize derivatives for f
 //  f->derivFDInit( 4 ); 

 //  // // Print the FD coefficients along the x, y, and z directions
 //  // const char *var[] = { "i", "j", "k" };
 //  // const Derivs::FiniteDiff::fd_t *fd_dd[] = { f->fd->fd_ddx, f->fd->fd_ddy, f->fd->fd_ddz }; 
  
 //  // for (int p=0; p<3; p++ ) {
 //  //   cout << "======================================================================" << endl;
 //  //   for (int i=0; i<f->dims[p]; i++ ) {
 //  //     cout << var[p] << " : " << i;
 //  //     for (int m=0; m<fd_dd[p][i].ssize; m++ ) {
 //  // 	cout << " " << fd_dd[p][i].coef[m];
 //  //     }
 //  //     cout << endl;
 //  //   }
 //  // }

 //  // const Derivs::FiniteDiff::fd_t *fd_dd2[] = { f->fd->fd_d2dx2, f->fd->fd_d2dy2, f->fd->fd_d2dz2 }; 
  
 //  // for (int p=0; p<3; p++ ) {
 //  //   cout << "======================================================================" << endl;
 //  //   for (int i=0; i<f->dims[p]; i++ ) {
 //  //     cout << var[p] << " : " << i;
 //  //     for (int m=0; m<fd_dd2[p][i].ssize; m++ ) {
 //  // 	cout << " " << fd_dd2[p][i].coef[m];
 //  //     }
 //  //     cout << endl;
 //  //   }
 //  // }

 //  MPIField *df = new MPIField( *f ); // Make a copy of f
 //  MPIField *df2 = new MPIField( *f ); // Make a copy of 

 //  // Assign the grid pointers; needed for finite differencing
 //  df->assignGrid( x, y, z );
 //  df2->assignGrid( x, y, z );
 //  // Initialize derivatives for f
 //  df->derivFDInit( 4 ); 
 //  df2->derivFDInit( 4 );
 //  // Take derivatives of f
 //  cout << "Taking derivatives" << endl;
 //  // df->ddx( *f );
 //  // df->ddy( *f );
 //  // df->ddz( *f );
 //  df->d2dx2( *f );
 //  df2->d2dy2( *f );
 //  // df->d2dz2( *f );
 //  // df->d2dxy( *f );
 //  // df->d2dxz( *f );
 //  // df->d2dyz( *f );

 //  df->add( *df, *df2 );

 //  // Set g to have same field data (does not copy entire class;
 //  // assumes they are initilized the same);
 //  (*g)=(*f);

 //  // Now take product of f and g and assign to f
 //  f->add(*f,*g);
 //  f->sub(*f,*g);
 //  f->mul(*f,*g);
 //  f->div(*f,*g); 

  if ( rank == 0 ) 
    {
      wtime = MPI_Wtime ( );

      printf ( "\n" );
      printf ( "HELLO_MPI - Master process:\n" );
      printf ( "  C/MPI version\n" );
      printf ( "  An MPI example program.\n" );
      printf ( "\n" );
      printf ( "  The number of processes is %d.\n", nproc );
      printf ( "\n" );
    }
  /*
  Every process prints a hello.
  */
  printf ( "  Process %d says 'Hello, world!'\n", rank );

  if ( rank == 0 )
    {
      wtime = MPI_Wtime ( ) - wtime;
      printf ( "  Elapsed wall clock time = %f seconds.\n", wtime );
    }

  MPI_Barrier( g->getMpiTopology()->comm );
  // Terminate MPI.
   MPI_Finalize ( );
  
  if ( rank == 0 )
    {
      printf ( "\n" );
      printf ( "HELLO_MPI - Master process:\n" );
      printf ( "  Normal end of execution: 'Goodbye, world!'\n" );
      printf ( "\n" );
    }

  return 0;
}
