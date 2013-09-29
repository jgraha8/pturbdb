#include <iostream>
//#include <stdlib.h>
#include <stdlib.h>
#include "esio/esio.h"
#include "mpi_topology.hpp"
#include "pturbdb_field.hpp"
#include "pfield_math.hpp"

#define DB_NZ 1536
#define DB_NY 512
#define DB_NX 2048

// #define FIELD_NZ 384
// #define FIELD_NY 256
// #define FIELD_NX 512

#define FIELD_NZ 64
#define FIELD_NY 64
#define FIELD_NX 64


using namespace std;
using namespace pturbdb;

int main(int argc, char *argv[]) {
	double wtime;

	MPI_Init(&argc, &argv);

	int db_dims[] = { DB_NZ, DB_NY, DB_NX };
	int db_field_offset[] = { FIELD_NZ, 0, FIELD_NX };
	int field_dims[] = { FIELD_NZ, FIELD_NY, FIELD_NX };
	int periodic[] = { 0, 0, 0 };

	PTurbDBField *u = new PTurbDBField("turbdb.conf", db_dims);

	u->PFieldInit(db_field_offset, field_dims, FIELD_DECOMP_PENCIL, periodic, 8);

	int rank = u->getMPITopology()->rank;
	int nproc = u->getMPITopology()->nproc;
	int *decomp_dims = u->getMPITopology()->dims;

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
		MPI_Barrier(u->getMPITopology()->comm);
	}

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

	if( u->getXLocal() == NULL ) {
		cout << "Grid not set as expected" << endl;
		int ierr=0;
		MPI_Abort( mpi_topology->comm, ierr);
	}

	MPI_Barrier(mpi_topology->comm);

	// Read the middle time
	double middle_time = 0.12345*0.5*( u->getDBTimeMax() + u->getDBTimeMin() );

	PTurbDBField *v = new PTurbDBField( *u, false ); // Do not copy u field data
	PTurbDBField *w = new PTurbDBField( *u, false ); // Do not copy u field data
	PTurbDBField *p = new PTurbDBField( *u, false ); // Do not copy u field data

	if( mpi_topology->rank == 0 ) cout << "Reading u from DB" << endl;
	// Read u from the DB
	//u->readDBField( 3*middle_time/2, "u" );
	*u = (double)10000.0L;
	if( mpi_topology->rank == 0 ) cout << "Reading v from DB" << endl;
	// Read u from the DB
	//v->readDBField( 3*middle_time/2, "v" );
	*v = (double)20.0;
	if( mpi_topology->rank == 0 ) cout << "Reading w from DB" << endl;
	// Read u from the DB
	//w->readDBField( 3*middle_time/2, "w" );
	*w = (double)-90.0;
	if( mpi_topology->rank == 0 ) cout << "Reading p from DB" << endl;
	// Read u from the DB
	p->readDBField( 3*middle_time/2, "p" );

	// // Read the operation grid
      	u->getXOperation( x_operation );
	u->getYOperation( y_operation );
	u->getZOperation( z_operation );

	// Create a velocity vector from the velocity components
	PFieldVector_t vel = PFieldVectorAssign( *u, *v, *w );

	// Velocity gradients
	if( mpi_topology->rank == 0 ) cout << "Computing u gradient" << endl;
	PFieldVector_t grad_u = PFieldGradient( *u );
	if( mpi_topology->rank == 0 ) cout << "Computing v gradient" << endl;
	PFieldVector_t grad_v = PFieldGradient( *v );
	if( mpi_topology->rank == 0 ) cout << "Computing w gradient" << endl;
	PFieldVector_t grad_w = PFieldGradient( *w );

	PFieldVector_t curl_vel = PFieldCurl( vel ) ;
	PFieldTensor_t grad_vel = PFieldTensorAssign( grad_u, grad_v, grad_w );
	
	// MPI_Barrier( mpi_topology->comm );
	// if( mpi_topology->rank == 0 ) cout << "Computing symmetric velocity gradient tensor" << endl;
	// PFieldTensor_t Sij = PFieldTensorSymmetric( grad_vel );

	// MPI_Barrier( mpi_topology->comm );
	// if( mpi_topology->rank == 0 ) cout << "Computing anti-symmetric velocity gradient tensor" << endl;
	// PFieldTensor_t Rij = PFieldTensorAntiSymmetric( grad_vel );

	// Create new fields 
	//PField *S = new PField( *u, false );
	//PField *R = new PField( *u, false );
	//PField *Q = new PField( *u, false );
	

	// MPI_Barrier( mpi_topology->comm );
	// if( mpi_topology->rank == 0 ) cout << "Computing scalar product of Sij" << endl;
	// PField *S = PFieldTensorDotDot( Sij, Sij );

	// MPI_Barrier( mpi_topology->comm );
	// if( mpi_topology->rank == 0 ) cout << "Computing scalar product of Rij" << endl;
	// PField *R = PFieldTensorDotDot( Rij, Rij );

	MPI_Barrier( mpi_topology->comm );
	if( mpi_topology->rank == 0 ) cout << "Computing Q invariant" << endl;
	PField *Q = PFieldTensorDotDot( grad_vel, grad_vel ); 
	*Q *= -0.5;

	// Temporary data buffer
	double *t_data = new double[u->getSizeOperation()];

	MPI_Barrier(mpi_topology->comm);
	if( mpi_topology->rank == 0 ) cout << "Writing output" << endl;

	// Output data
	esio_handle h = esio_handle_initialize(mpi_topology->comm);

	// Open the database file
	esio_file_create(h, "/datascope/tdbchannel/analysis/q.h5", 1); // Open read-only

	esio_field_establish(h, field_dims[0], offset_local[0]+offset_operation[0], dims_operation[0],
			        field_dims[1], offset_local[1]+offset_operation[1], dims_operation[1],
			        field_dims[2], offset_local[2]+offset_operation[2], dims_operation[2]);

	if( mpi_topology->rank == 0 ) cout << "Writing grid" << endl;
	
	double *x_full = new double[u->getSizeOperation()];
	double *y_full = new double[u->getSizeOperation()];
	double *z_full = new double[u->getSizeOperation()];

	long index=0;
	for(int i=0; i<dims_operation[0]; i++) {
		for (int j=0; j<dims_operation[1]; j++) {
			for(int k=0; k<dims_operation[2]; k++) {
				x_full[index] = x_operation[i];
				y_full[index] = y_operation[j];
				z_full[index] = z_operation[k];
				index++;
			}
		}
	}
	esio_field_write_double(h, grid_field_names[0], x_full, 0, 0, 0, "z");
	esio_field_write_double(h, grid_field_names[1], y_full, 0, 0, 0, "y");
	esio_field_write_double(h, grid_field_names[2], z_full, 0, 0, 0, "x");

	if( mpi_topology->rank == 0 ) cout << "Writing fields" << endl;

	u->getDataOperation(t_data);
	esio_field_write_double(h, "u", t_data, 0, 0, 0, "u");

	v->getDataOperation(t_data);
	esio_field_write_double(h, "v", t_data, 0, 0, 0, "v");

	w->getDataOperation(t_data);
	esio_field_write_double(h, "w", t_data, 0, 0, 0, "w");

	p->getDataOperation(t_data);
	esio_field_write_double(h, "p", t_data, 0, 0, 0, "p");

	curl_vel[0]->getDataOperation(t_data);
	esio_field_write_double(h, "omega_x", t_data, 0, 0, 0, "omega_x");
	curl_vel[1]->getDataOperation(t_data);
	esio_field_write_double(h, "omega_y", t_data, 0, 0, 0, "omega_y");
	curl_vel[2]->getDataOperation(t_data);
	esio_field_write_double(h, "omega_z", t_data, 0, 0, 0, "omega_z");

	Q->getDataOperation(t_data);
	esio_field_write_double(h, "Q", t_data, 0, 0, 0, "Q");

	// dudx->getDataOperation(t_data);
	// esio_field_write_double(h, "dudx", t_data, 0, 0, 0, "dudx");

	// dudy->getDataOperation(t_data);
	// esio_field_write_double(h, "dudy", t_data, 0, 0, 0, "dudy");

	// dudz->getDataOperation(t_data);
	// esio_field_write_double(h, "dudz", t_data, 0, 0, 0, "dudz");

	// dvdx->getDataOperation(t_data);
	// esio_field_write_double(h, "dvdx", t_data, 0, 0, 0, "dvdx");

	// dvdy->getDataOperation(t_data);
	// esio_field_write_double(h, "dvdy", t_data, 0, 0, 0, "dvdy");

	// dvdz->getDataOperation(t_data);
	// esio_field_write_double(h, "dvdz", t_data, 0, 0, 0, "dvdz");

	// dwdx->getDataOperation(t_data);
	// esio_field_write_double(h, "dwdx", t_data, 0, 0, 0, "dwdx");

	// dwdy->getDataOperation(t_data);
	// esio_field_write_double(h, "dwdy", t_data, 0, 0, 0, "dwdy");

	// dwdz->getDataOperation(t_data);
	// esio_field_write_double(h, "dwdz", t_data, 0, 0, 0, "dwdz");
	
	esio_file_close(h);
	esio_handle_finalize(h);


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

	MPI_Barrier(mpi_topology->comm);

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
