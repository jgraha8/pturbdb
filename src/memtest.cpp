#include <iostream>
#include <sstream>
//#include <stdlib.h>
#include <stdlib.h>
#include "esio/esio.h"
#include "mpi_topology.hpp"
#include "pturbdb_field.hpp"
#include "pfield_math.hpp"

#define DB_NZ 1536
#define DB_NY 512
#define DB_NX 2048

#define DB_DT 0.0065
// #define FIELD_NZ 384
// #define FIELD_NY 256
// #define FIELD_NX 512

#define FIELD_NZ 32
#define FIELD_NY 32
#define FIELD_NX 32


using namespace std;
using namespace pturbdb;

int main(int argc, char *argv[]) {
	double wtime;

	MPI_Init(&argc, &argv);

	int db_dims[] = { DB_NZ, DB_NY, DB_NX };
	int db_field_offset[] = { FIELD_NZ, (DB_NY-FIELD_NY)/2, FIELD_NX };
 	int field_dims[] = { FIELD_NZ, FIELD_NY, FIELD_NX };
	int periodic[] = { 0, 0, 0 };

	PTurbDBField *u = new PTurbDBField("turbdb.conf", db_dims);

	u->PFieldInit(db_field_offset, field_dims, FIELD_DECOMP_PENCIL, periodic, 4);

	int rank = u->getMPITopology()->rank;
	int nproc = u->getMPITopology()->nproc;
	int *decomp_dims = u->getMPITopology()->dims;

	// Read the grid
	const char *grid_field_names[] = { "z", "y", "x" };

	int *dims_local = u->getDimsLocal();
	int *dims_operation = u->getDimsOperation();
	int *offset_local = u->getOffsetLocal();
	int *offset_operation = u->getOffsetOperation();
	MPITopology_t *mpi_topology = u->getMPITopology();

	double *x = new double[dims_local[0]];
	double *y = new double[dims_local[1]];
	double *z = new double[dims_local[2]];

	double *x_operation = new double[dims_operation[0]];
	double *y_operation = new double[dims_operation[1]];
	double *z_operation = new double[dims_operation[2]];

	// This sets the grid
	u->readDBGridLocal(grid_field_names, x, y, z);

	MPI_Barrier(mpi_topology->comm);

	// // Read the operation grid
      	u->getXOperation( x_operation );
	u->getYOperation( y_operation );
	u->getZOperation( z_operation );

	// Set the starting time and the time step
	const double start_time = u->getDBTimeMax() - DB_DT * 43;
	static const double dt = DB_DT * 2;

	//PTurbDBField *v = new PTurbDBField( *u, false ); // Do not copy u field data
	//	PTurbDBField *w = new PTurbDBField( *u, false ); // Do not copy u field data
	//	PTurbDBField *p = new PTurbDBField( *u, false ); // Do not copy u field data

	if( mpi_topology->rank == 0 ) cout << "Reading u from DB" << endl;
	u->readDBField( start_time, "u" );

	// if( mpi_topology->rank == 0 ) cout << "Reading v from DB" << endl;
	// v->readDBField( time, "v" );

	// if( mpi_topology->rank == 0 ) cout << "Reading w from DB" << endl;
	// w->readDBField( time, "w" );

	// if( mpi_topology->rank == 0 ) cout << "Reading p from DB" << endl;
	// p->readDBField( time, "p" );

	// for (int n=0; n<-1; n++) {

	// 	double time = start_time - dt * n;

	// 	// Create a velocity vector from the velocity components
	// 	PFieldVector_t vel = PFieldVectorAssign( u, v, w );

	// 	// Velocity gradients
	// 	if( mpi_topology->rank == 0 ) cout << "Computing u gradient" << endl;
	// 	PFieldVector_t grad_u = PFieldGradient( *u );
	// 	if( mpi_topology->rank == 0 ) cout << "Computing v gradient" << endl;
	// 	PFieldVector_t grad_v = PFieldGradient( *v );
	// 	if( mpi_topology->rank == 0 ) cout << "Computing w gradient" << endl;
	// 	PFieldVector_t grad_w = PFieldGradient( *w );

	// 	PFieldTensor_t Aij = PFieldTensorAssign( grad_u, grad_v, grad_w );
	// 	PFieldTensor_t Aji = PFieldTensorNew( Aij );
	// 	// Clear vectors and tensors created with *Assign
	// 	vel.clear();
	// 	Aij.clear();

	// 	PFieldVectorDelete( grad_u );
	// 	PFieldVectorDelete( grad_v );
	// 	PFieldVectorDelete( grad_w );

	// 	PFieldTensorDelete( Aji );
 
	// }

	delete u;
	//	delete v;
	//	delete w;
	//	delete p;

	delete [] x;
	delete [] y;
	delete [] z;
	delete [] x_operation;
	delete [] y_operation;
	delete [] z_operation;

	// Terminate MPI.
	MPI_Finalize();

	return 0;
}
