#include <iostream>
//#include <stdlib.h>
#include <stdlib.h>

#include "turbdb_field.hpp"

#define DB_NZ 1536
#define DB_NY 512
#define DB_NX 2048

#define FIELD_NZ 256
#define FIELD_NY 256
#define FIELD_NX 256

//#define FIELD_NZ 128
//#define FIELD_NY 64
//#define FIELD_NX 128

using namespace std;
using namespace pturb_fields;

int main(int argc, char *argv[]) {
	double wtime;

	MPI_Init(&argc, &argv);

	int db_dims[] = { DB_NZ, DB_NY, DB_NX };
	int db_field_offset[] = { 0, 0, 0 };
	int field_dims[] = { FIELD_NZ, FIELD_NY, FIELD_NX };
	int periodic[] = { 1, 1, 1 };

	TurbDBField *u = new TurbDBField("turbdb.conf", db_dims);

	u->dbFieldInit(db_field_offset, field_dims, FIELD_DECOMP_PENCIL, periodic, 4);

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

	// This sets the grid
	u->readDBGridLocal(grid_field_names, x, y, z);

	if( u->getXLocal() == NULL ) {
		cout << "Grid not set as expected" << endl;
		int ierr=0;
		MPI_Abort( u->getMpiTopology()->comm, ierr);
	}
		

	MPI_Barrier(u->getMpiTopology()->comm);
	// // Print the grid to stdout
	// for (int n = 0; n < nproc; n++) {
	// 	if (n == rank) {
	// 		for (long i = 0; i < dims_local[2]; i++) {
	// 			cout << "rank, i, z[i] " << rank << " " << i << " " << z[i]
	// 					<< endl;
	// 		}
	// 	}
	// 	MPI_Barrier(u->getMpiTopology()->comm);
	// }

	// Read the middle time
	double middle_time = 0.5*( u->getDBTimeMax() + u->getDBTimeMin() );

	TurbDBField *v = new TurbDBField( *u, false ); // Do not copy u field data
	TurbDBField *w = new TurbDBField( *u, false ); // Do not copy u field data

	// Read u from the DB
	u->readDBField( middle_time, "u" );
	// Read u from the DB
	v->readDBField( middle_time, "v" );
	// Read u from the DB
	w->readDBField( middle_time, "w" );

	Field *dudx = new Field( field_dims, FIELD_DECOMP_SLAB, periodic, 4 );
	// Set the grid and initialize finite differences
	dudx->setGridLocal( x, y, z );
	dudx->finiteDiffInit();

	Field *dudy = new Field( *dudx );
	Field *dudz = new Field( *dudx );

	Field *dvdx = new Field( *dudx );
	Field *dvdy = new Field( *dudx );
	Field *dvdz = new Field( *dudx );

	Field *dwdx = new Field( *dudx );
	Field *dwdy = new Field( *dudx );
	Field *dwdz = new Field( *dudx );

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

	Field *Q = new Field( *dudx, false );
	Field *S11 = new Field( *dudx );
	Field *S12 = new Field( *dudy );
	Field *S13 = new Field( *dudz );
	Field *S22 = new Field( *dvdy );
	Field *S23 = new Field( *dvdz );
	Field *S33 = new Field( *dwdz );

	Field *O12 = new Field( *dudy );
	Field *O13 = new Field( *dudz );
	Field *O23 = new Field( *dvdz );

	(*S12 += (*dvdx))*=0.5;
	(*S13 += (*dwdx))*=0.5;
	(*S23 += (*dwdy))*=0.5;

	(*O12 -= (*dvdx))*=0.5;
	(*O13 -= (*dwdx))*=0.5;
	(*O23 -= (*dwdy))*=0.5;

	Field *TrS = new Field( *S11, false );
	Field *TrO = new Field( *S11, false );

	// Here using Q as a buffer. The multiplication takes place and is stored in Q. We then add it to TrS.
	TrO->mul( *S12, *S12);
	*TrO += Q->mul( *S13, *S13 );
	*TrO += Q->mul( *S23, *S23 );
	*TrO *= -2.0; // Multiply the result by -2.0

	TrS->mul( *S11, *S11 );
	*TrS += Q->mul( *S22, *S22 );
	*TrS += Q->mul( *S33, *S33 );
	*TrS -= *TrO; // Subtract the trace of the antisymmetric component; they have the same parts

	// Q then only retains 1/2 of the diagonal of the symmetric component
	Q->sub( *TrO, *TrS ) *= 0.5;

	

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
