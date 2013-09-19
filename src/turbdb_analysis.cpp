#include <iostream>
//#include <stdlib.h>
#include <stdlib.h>

#include "turbdb_field.hpp"

#define DB_NZ 1536
#define DB_NY 512
#define DB_NX 2048

//#define FIELD_NZ 1536
//#define FIELD_NY 64
//#define FIELD_NX 2048

#define FIELD_NZ 128
#define FIELD_NY 64
#define FIELD_NX 128

using namespace std;
using namespace pturb_fields;

int main(int argc, char *argv[]) {
	double wtime;

	MPI_Init(&argc, &argv);

	int ndim = 3;
	int db_dims[] = { DB_NZ, DB_NY, DB_NX };
	int db_field_offset[] = { 0, 0, 0 };
	int field_dims[] = { FIELD_NZ, FIELD_NY, FIELD_NX };
	int periodic[] = { 0, 0, 0 };

	TurbDBField *u = new TurbDBField("turbdb.conf", db_dims);

	u->dbFieldInit(db_field_offset, field_dims, FIELD_DECOMP_PENCIL, periodic,
			4);

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

	double *x = new double[dims_local[0]];
	double *y = new double[dims_local[1]];
	double *z = new double[dims_local[2]];

	u->readDBGridLocal(grid_field_names, x, y, z);

	MPI_Barrier(u->getMpiTopology()->comm);
	// Print the grid to stdout
	for (int n = 0; n < nproc; n++) {
		if (n == rank) {
			for (long i = 0; i < dims_local[2]; i++) {
				cout << "rank, i, z[i] " << rank << " " << i << " " << z[i]
						<< endl;
			}
		}
		MPI_Barrier(u->getMpiTopology()->comm);
	}

	// Read the middle time
	double middle_time = 0.5*( u->getDBTimeMax() + u->getDBTimeMin() );

	TurbDBField *v = new TurbDBField( *u );
	TurbDBField *w = new TurbDBField( *u );

	v->dbFieldInit(db_field_offset, field_dims, FIELD_DECOMP_PENCIL, periodic, 4);
	w->dbFieldInit(db_field_offset, field_dims, FIELD_DECOMP_PENCIL, periodic, 4);

	// Assign the grid pointers; needed for finite differencing
	// u was set on the readDBGridLocal call
	v->setGridLocal( x, y, z );
	w->setGridLocal( x, y, z );

	//  // Initialize derivatives for f
	u->finiteDiffInit();
	v->finiteDiffInit();
	w->finiteDiffInit();

	// Read u from the DB
	u->readDBField( middle_time, "u" );
	// Read u from the DB
	v->readDBField( middle_time, "v" );
	// Read u from the DB
	w->readDBField( middle_time, "w" );

	Field *dudx = new Field( field_dims, FIELD_DECOMP_PENCIL, periodic, 4 );
	Field *dudy = new Field( *dudx );
	Field *dudz = new Field( *dudz );

	Field *dvdx = new Field( *dudx );
	Field *dvdy = new Field( *dudx );
	Field *dvdz = new Field( *dudz );

	Field *dwdx = new Field( *dudx );
	Field *dwdy = new Field( *dudz );
	Field *dwdz = new Field( *dudx );

	dudx->setGridLocal( x, y, z );
	dudy->setGridLocal( x, y, z );
	dudz->setGridLocal( x, y, z );

	dvdx->setGridLocal( x, y, z );
	dvdy->setGridLocal( x, y, z );
	dvdz->setGridLocal( x, y, z );

	dwdx->setGridLocal( x, y, z );
	dwdy->setGridLocal( x, y, z );
	dwdz->setGridLocal( x, y, z );

	dudx->finiteDiffInit();
	dudy->finiteDiffInit();
	dudz->finiteDiffInit();

	dvdx->finiteDiffInit();
	dvdy->finiteDiffInit();
	dvdz->finiteDiffInit();

	dwdx->finiteDiffInit();
	dwdy->finiteDiffInit();
	dwdz->finiteDiffInit();

	dudx->ddx( *u );
	dudy->ddy( *u );
	dudz->ddz( *u );

	dvdx->ddx( *v );
	dvdy->ddy( *v );
	dvdz->ddz( *v );

	dwdx->ddx( *w );
	dwdy->ddy( *w );
	dwdz->ddz( *w );

	//Field *Q = new Field( *dudx );
	

	// Compute the Q-criterion
	
	MPI_Barrier(test->getMpiTopology()->comm);

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
