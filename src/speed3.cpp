#include <iostream>
//#include <stdlib.h>
#include <stdlib.h>

#include "field.hpp"

#define NX 32
#define NY 32
#define NZ 32

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
  int periodic[] = {0, 0, 0};

  Field *f = new Field( dims, FIELD_DECOMP_PENCIL, periodic, 2 );
  Field *g = new Field( *f ); // Make a copy of f

  int rank = f->getMpiTopology()->rank;
  int nproc = f->getMpiTopology()->nproc;
  int *decomp_dims = f->getMpiTopology()->dims;

  for (int n=0; n<nproc; n++) {
    if( n == rank ) {
      cout << "Number of global data elements : " << g->getSize() << endl;
      cout << "Number of local data elements : " << g->getSizeLocal() << endl;
      cout << "Number of operation data elements : " << g->getSizeOperation() << endl;
      cout << "MPI Topology : " << decomp_dims[0] << "x" << decomp_dims[1] << "x" << decomp_dims[2] << endl;
    }
    MPI_Barrier( g->getMpiTopology()->comm );
  }

  int *dims_local = f->getDimsLocal();

  long index=0;
  for (int i=0; i<dims_local[0]; i++) {
    for (int j=0; j<dims_local[1]; j++) {
	for (int k=0; k<dims_local[2]; k++) {
	    f->data_local[index++] = (double)i*j + k;
	}
    }
  }
  
  double *x = new double[*dims_local];
  double *y = new double[*(dims_local+1)];
  double *z = new double[*(dims_local+2)];

  for (int i=0; i<dims_local[0]; i++ ) x[i] = (double)i;
  for (int j=0; j<dims_local[1]; j++ ) y[j] = (double)j;
  for (int k=0; k<dims_local[2]; k++ ) z[k] = (double)k;

  // Assign the grid pointers; needed for finite differencing
  f->setGridLocal( x, y, z );

  // Initialize derivatives for f
  f->finiteDiffInit(); 
       
  Field *df = new Field( *f ); // Make a copy of f
  //  MPIField *df2 = new MPIField( *f ); // Make a copy of 

  //  // Assign the grid pointers; needed for finite differencing
  df->setGridLocal( x, y, z );
  //  df2->assignGrid( x, y, z );
  //  // Initialize derivatives for f
  df->finiteDiffInit(); 
  //  df2->derivFDInit( 4 );
  // Take derivatives of f
  cout << "Taking derivatives" << endl;
  df->ddx( *f );
  *f+=*f;
  df->ddy( *f );
  *f-=*f;
  df->ddz( *f );
  df->d2dx2( *f );
  if( f->getSynchronized() == false ) cout << "something bad" << endl;
  df->d2dy2( *f );
  df->d2dz2( *f );
  df->d2dxy( *f );
  df->d2dxz( *f );
  df->d2dyz( *f );

  // Test synchronization
  Field *test = new Field( *f );

  dims_local = test->getDimsLocal();
  int *dims_operation = test->getDimsOperation();

  index=0;
  for (int i=0; i<dims_local[0]; i++) {
    for (int j=0; j<dims_local[1]; j++) {
      for (int k=0; k<dims_local[2]; k++) {
	if( test->getMpiTopology()->coords[0] == 0 && 
	    test->getMpiTopology()->coords[1] == 0 ) {
	  test->data_local[index++] = 0.0;
	} else {
	  test->data_local[index++] = 1.0;
	}
      }
    }
  }

  //  test->synchronize();
  df->ddx( *test );

  double sum_test=0.0;
  if( test->getMpiTopology()->coords[0] == 0 && 
      test->getMpiTopology()->coords[1] == 0 ) {

    cout << "center rank : " << test->getMpiTopology()->rank << endl;
    cout << "dims_local : " << dims_local[0] << " " << dims_local[1] << " " << dims_local[2] << endl;
    cout << "dims_operation : " << dims_operation[0] << " " << dims_operation[1] << " " << dims_operation[2] << endl;

    index=0;
    sum_test=0.0;
    for (int i=0; i<dims_local[0]; i++) {
      for (int j=0; j<dims_local[1]; j++) {
	for (int k=0; k<dims_local[2]; k++) {
	  sum_test += test->data_local[index++];
	}
      }
    }

    // sum test should equal the number of rind points
    cout << "sum_test, # rind points : " << (long)sum_test << " ";
    cout << test->getSizeRind(0,-1) + test->getSizeRind(0,1) + test->getSizeRind(1,-1) + test->getSizeRind(1,1) << endl;

  }

  MPI_Barrier(test->getMpiTopology()->comm);

  
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
