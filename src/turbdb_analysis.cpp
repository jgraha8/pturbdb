#include <iostream>
//#include <stdlib.h>
#include <stdlib.h>

#include "pturbdb_field.hpp"

#define DB_NZ 1536
#define DB_NY 512
#define DB_NX 2048

#define FIELD_NZ 64
#define FIELD_NY 64
#define FIELD_NX 64

//#define FIELD_NZ 128
//#define FIELD_NY 64
//#define FIELD_NX 128

using namespace std;
using namespace pturbdb;

int main(int argc, char *argv[]) {
	double wtime;

	MPI_Init(&argc, &argv);

	int db_dims[] = { DB_NZ, DB_NY, DB_NX };
	int db_field_offset[] = { 345, 127, 789 };
	int field_dims[] = { FIELD_NZ, FIELD_NY, FIELD_NX };
	int periodic[] = { 0, 0, 0 };

	PTurbDBField *u = new PTurbDBField("turbdb.conf", db_dims);

	u->PFieldInit(db_field_offset, field_dims, FIELD_DECOMP_PENCIL, periodic, 4);

	int rank = u->getMpiTopology()->rank;
	int nproc = u->getMpiTopology()->nproc;
	int *decomp_dims = u->getMpiTopology()->dims;

	for (int n = 0; n < nproc; n++) {
		if (n == rank) {
			cout << "Number of global data elements : " << u->getSize() << endl;
			cout << "Number of local data elements : " << u->getSizeLocal()
					<< endl;
			cout << "Number of operation data elements : "
					<< u->getSizeOperation() << endl;
			cout << "MPI Topology : " << decomp_dims[0] << "x" << decomp_dims[1]
					<< "x" << decomp_dims[2] << endl;
			cout << endl;
			cout << "DB time range : " << u->getDBTimeMin() << "-"
					<< u->getDBTimeMax() << endl;

		}
		MPI_Barrier(u->getMpiTopology()->comm);
	}

	// Read the grid
	const char *grid_field_names[] = { "z", "y", "x" };

	int *dims_local = u->getDimsLocal();
	int *dims_operation = u->getDimsOperation();

	double *x = new double[dims_local[0]];
	double *y = new double[dims_local[1]];
	double *z = new double[dims_local[2]];

	double *x_operation = new double[dims_operation[0]];
	double *y_operation = new double[dims_operation[1]];
	double *z_operation = new double[dims_operation[2]];

	// This sets the grid
	u->readDBGridLocal(grid_field_names, x, y, z);

	if( u->getXLocal() == NULL ) {
		cout << "Grid not set as expected" << endl;
		int ierr=0;
		MPI_Abort( u->getMpiTopology()->comm, ierr);
	}

	// Read the operation grid
      	u->getXOperation( x_operation );
	u->getYOperation( y_operation );
	u->getZOperation( z_operation );

		

	MPI_Barrier(u->getMpiTopology()->comm);

	// Read the middle time
	double middle_time = 0.5*( u->getDBTimeMax() + u->getDBTimeMin() );

	PTurbDBField *v = new PTurbDBField( *u, false ); // Do not copy u field data
	PTurbDBField *w = new PTurbDBField( *u, false ); // Do not copy u field data

	// Read u from the DB
	u->readDBField( middle_time/2, "u" );
	// Read u from the DB
	v->readDBField( middle_time/2, "v" );
	// Read u from the DB
	w->readDBField( middle_time/2, "w" );

	PField *dudx = new PField( field_dims, FIELD_DECOMP_SLAB, periodic, 4 );
	// Set the grid and initialize finite differences
	dudx->setGridLocal( x, y, z );
	dudx->finiteDiffInit();

	PField *dudy = new PField( *dudx );
	PField *dudz = new PField( *dudx );

	PField *dvdx = new PField( *dudx );
	PField *dvdy = new PField( *dudx );
	PField *dvdz = new PField( *dudx );

	PField *dwdx = new PField( *dudx );
	PField *dwdy = new PField( *dudx );
	PField *dwdz = new PField( *dudx );

	if( u->getMpiTopology()->rank == 0 ) cout << "Computing u derivatives" << endl;

	dudx->ddx( *u );
	dudy->ddy( *u );
	dudz->ddz( *u );

	if( u->getMpiTopology()->rank == 0 ) cout << "Computing v derivatives" << endl;

	dvdx->ddx( *v );
	dvdy->ddy( *v );
	dvdz->ddz( *v );

	if( u->getMpiTopology()->rank == 0 ) cout << "Computing w derivatives" << endl;

	dwdx->ddx( *w );
	dwdy->ddy( *w );
	dwdz->ddz( *w );

	PField *Q = new PField( *dudx, false );
	PField *S11 = new PField( *dudx );
	PField *S12 = new PField( *dudy );
	PField *S13 = new PField( *dudz );
	PField *S22 = new PField( *dvdy );
	PField *S23 = new PField( *dvdz );
	PField *S33 = new PField( *dwdz );

	PField *O12 = new PField( *dudy );
	PField *O13 = new PField( *dudz );
	PField *O23 = new PField( *dvdz );

	(*S12 += *dvdx)*=0.5;
	(*S13 += *dwdx)*=0.5;
	(*S23 += *dwdy)*=0.5;

	(*O12 -= *dvdx)*=0.5;
	(*O13 -= *dwdx)*=0.5;
	(*O23 -= *dwdy)*=0.5;

	PField *TrS = new PField( *S11, false );
	PField *TrO = new PField( *S11, false );

	// Here using Q as a buffer. The multiplication takes place and is stored in Q. We then add it to TrS.
	// TrO->mul( *S12, *S12);
	// *TrO += Q->mul( *S13, *S13 );
	// *TrO += Q->mul( *S23, *S23 );
	// *TrO *= -2.0; // Multiply the result by -2.0

	*TrO = Q->add( Q->add( Q->mul( *S12, 
				       *S12 ),
			       Q->mul( *S13, 
				       *S13 ) ), 
		       Q->mul( *S23, 
			       *S23 ) ) *= -2.0;

	TrS->mul( *S11, *S11 );
	*TrS += Q->mul( *S22, *S22 );
	*TrS += Q->mul( *S33, *S33 );
	*TrS -= *TrO; // Subtract the trace of the antisymmetric component; they have the same parts

	// Q then only retains 1/2 of the diagonal of the symmetric component
	Q->sub( *TrO, *TrS ) *= 0.5;

	double *Q_data = new double[Q->getSizeOperation()];

	Q->getDataOperation( Q_data );
	// Now set the data
	*Q = Q_data;

	// Compute the Q-criterion
	
	MPI_Barrier(u->getMpiTopology()->comm);

	// //  df->add( *df, *df2 );

	// //  // Set g to have same field data (does not copy entire class;
	// //  // assumes they are initilized the same);
	// //  (*g)=(*f);

	// //  // Now take product of f and g and assign to f
	// //  f->add(*f,*g);
	// //  f->sub(*f,*g);
	// //  f->mul(*f,*g);
	// //  f->div(*f,*g);
	
	if (rank == 0) {
		wtime = MPI_Wtime();

		printf("\n");
		printf("HELLO_MPI - Master process:\n");
		printf("  C/MPI version\n");
		printf("  An MPI example program.\n");
		printf("\n");
		printf("  The number of processes is %d.\n", nproc);
		printf("\n");
	}
	/*
	 Every process prints a hello.
	 */
	printf("  Process %d says 'Hello, world!'\n", rank);

	if (rank == 0) {
		wtime = MPI_Wtime() - wtime;
		printf("  Elapsed wall clock time = %f seconds.\n", wtime);
	}

	MPI_Barrier(u->getMpiTopology()->comm);

	// Terminate MPI.
	MPI_Finalize();

	if (rank == 0) {
		printf("\n");
		printf("HELLO_MPI - Master process:\n");
		printf("  Normal end of execution: 'Goodbye, world!'\n");
		printf("\n");
	}

	return 0;
}
