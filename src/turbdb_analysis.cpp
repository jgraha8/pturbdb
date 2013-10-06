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

#define FIELD_NZ 128
#define FIELD_NY 128
#define FIELD_NX 128


using namespace std;
using namespace pturbdb;

int main(int argc, char *argv[]) {
	double wtime;

	MPI_Init(&argc, &argv);

	const int db_dims[] = { DB_NZ, DB_NY, DB_NX };
	//	int db_field_offset[] = { FIELD_NZ, (DB_NY-FIELD_NY)/2, FIELD_NX };
	const int db_field_offset[] = { FIELD_NZ, 0, FIELD_NX };
	const int field_dims[] = { FIELD_NZ, FIELD_NY, FIELD_NX };
	const int periodic[] = { 0, 0, 0 };

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

	// // Read the operation grid
      	u->getXOperation( x_operation );
	u->getYOperation( y_operation );
	u->getZOperation( z_operation );

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
	
	// Set the starting time and the time step
	const double start_time = u->getDBTimeMax()/2;
	static const double dt = DB_DT/2;

	PTurbDBField *v = new PTurbDBField( *u, false ); // Do not copy u field data
	PTurbDBField *w = new PTurbDBField( *u, false ); // Do not copy u field data
	PTurbDBField *p = new PTurbDBField( *u, false ); // Do not copy u field data

	PFieldVector_t vel    = PFieldVectorNew( *u );
	PFieldVector_t grad_u = PFieldVectorNew( *u );
	PFieldVector_t grad_v = PFieldVectorNew( *u );
	PFieldVector_t grad_w = PFieldVectorNew( *u );
	PFieldVector_t vorticity = PFieldVectorNew( *u );

	PFieldTensor_t Aij = PFieldTensorNew( *u );
	PFieldTensor_t Sij = PFieldTensorNew( *u );
	PFieldTensor_t Tij = PFieldTensorNew( *u );
	PFieldTensor_t Tji = PFieldTensorNew( *u );
	
	PField *S2 = new PField( *u, false );
	PField *T2 = new PField( *u, false );
	PField *Q  = new PField( *u, false );
	PField *Q1 = new PField( *u, false );

	// Temporary data buffer
	double *t_data = new double[u->getSizeOperation()];

	static const int nsteps = 300;
	for (int n=0; n<nsteps; n++) {

		if( mpi_topology->rank == 0 ) cout << "Percent complete " << n*100.0/nsteps << endl;

		double time = start_time - dt * n;

		if( mpi_topology->rank == 0 ) cout << "Reading u from DB" << endl;
		u->readDBField( time, "u" );

		if( mpi_topology->rank == 0 ) cout << "Reading v from DB" << endl;
		v->readDBField( time, "v" );

		if( mpi_topology->rank == 0 ) cout << "Reading w from DB" << endl;
		w->readDBField( time, "w" );

		if( mpi_topology->rank == 0 ) cout << "Reading p from DB" << endl;
		p->readDBField( time, "p" );

		// Assign pointers to the velocity vector
		PFieldVectorAssign( vel, u, v, w );

		// Velocity gradients
		if( mpi_topology->rank == 0 ) cout << "Computing u gradient" << endl;
		PFieldVectorGradient( grad_u, *u );
		if( mpi_topology->rank == 0 ) cout << "Computing v gradient" << endl;
		PFieldVectorGradient( grad_v, *v );
		if( mpi_topology->rank == 0 ) cout << "Computing w gradient" << endl;
		PFieldVectorGradient( grad_w, *w );

		PFieldVectorCurl( vorticity, vel ) ;
		PFieldTensorAssign( Aij, grad_u, grad_v, grad_w );

		MPI_Barrier( mpi_topology->comm );
		if( mpi_topology->rank == 0 ) cout << "Computing symmetric velocity gradient tensor" << endl;
		PFieldTensorSymmetric( Sij, Aij );

		MPI_Barrier( mpi_topology->comm );
		if( mpi_topology->rank == 0 ) cout << "Computing anti-symmetric velocity gradient tensor" << endl;
		PFieldTensorAntiSymmetric( Tij, Aij );
		
		MPI_Barrier( mpi_topology->comm );
		if( mpi_topology->rank == 0 ) cout << "Computing transpose of anti-symmetric velocity gradient tensor" << endl;
		PFieldTensorTranspose( Tji, Tij );
		
		MPI_Barrier( mpi_topology->comm );
		if( mpi_topology->rank == 0 ) cout << "Computing scalar product of Sij" << endl;
		PFieldTensorDotDot( *S2, Sij, Sij );

		MPI_Barrier( mpi_topology->comm );
		if( mpi_topology->rank == 0 ) cout << "Computing scalar product of Rij" << endl;
		PFieldTensorDotDot( *T2, Tij, Tji );

		MPI_Barrier( mpi_topology->comm );

		if( mpi_topology->rank == 0 ) cout << "Computing Q invariant" << endl;
		//PField *Q = PFieldTensorDotDot( Aij, Aij ); 
		// *Q *= -(double)0.5L; // Apply the -1/2 factor

		Q->sub( *T2, *S2 ) *= 0.5;
		
		PFieldTensorDotDot( *Q1, Aij, Aij ) *= -(double)0.5; // Apply the -1/2 factor

		MPI_Barrier(mpi_topology->comm);
		if( mpi_topology->rank == 0 ) cout << "Writing output" << endl;
		
		// Output data
		esio_handle h = esio_handle_initialize(mpi_topology->comm);

		//std::stringstream *h5file = new std::stringstream;
		// *h5file << "/datascope/tdbchannel/analysis/q-wall-time-";
		// h5file->precision(8);
		// *h5file << time << ".h5";
		
		std::string *h5file = new std::string("/datascope/tdbchannel/analysis/q-wall-time-");
		char time_buffer[10];
		sprintf(time_buffer, "%9.6f", time);
		char *s = time_buffer;
		while( *s == ' ' ) *(s++)='0';


		*h5file = *h5file + time_buffer + ".h5";

		// Open the database file
		esio_file_create(h, h5file->c_str(), 1);

		delete h5file;
		
		esio_field_establish(h, field_dims[0], offset_local[0]+offset_operation[0], dims_operation[0],
				     field_dims[1], offset_local[1]+offset_operation[1], dims_operation[1],
				     field_dims[2], offset_local[2]+offset_operation[2], dims_operation[2]);

		if( mpi_topology->rank == 0 ) cout << "Writing grid" << endl;
	
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

		vorticity[0]->getDataOperation(t_data);
		esio_field_write_double(h, "omega_x", t_data, 0, 0, 0, "omega_x");
		vorticity[1]->getDataOperation(t_data);
		esio_field_write_double(h, "omega_y", t_data, 0, 0, 0, "omega_y");
		vorticity[2]->getDataOperation(t_data);
		esio_field_write_double(h, "omega_z", t_data, 0, 0, 0, "omega_z");

		Q->getDataOperation(t_data);
		esio_field_write_double(h, "Q", t_data, 0, 0, 0, "Q");

		Q1->getDataOperation(t_data);
		esio_field_write_double(h, "Q1", t_data, 0, 0, 0, "Q1");

		// Aij[0][0]->getDataOperation(t_data);
		// esio_field_write_double(h, "dudx", t_data, 0, 0, 0, "dudx");

		// Aij[0][1]->getDataOperation(t_data);
		// esio_field_write_double(h, "dudy", t_data, 0, 0, 0, "dudy");

		// Aij[0][2]->getDataOperation(t_data);
		// esio_field_write_double(h, "dudz", t_data, 0, 0, 0, "dudz");

		// Aij[1][0]->getDataOperation(t_data);
		// esio_field_write_double(h, "dvdx", t_data, 0, 0, 0, "dvdx");

		// Aij[1][1]->getDataOperation(t_data);
		// esio_field_write_double(h, "dvdy", t_data, 0, 0, 0, "dvdy");

		// Aij[1][2]->getDataOperation(t_data);
		// esio_field_write_double(h, "dvdz", t_data, 0, 0, 0, "dvdz");

		// Aij[2][0]->getDataOperation(t_data);
		// esio_field_write_double(h, "dwdx", t_data, 0, 0, 0, "dwdx");

		// Aij[2][1]->getDataOperation(t_data);
		// esio_field_write_double(h, "dwdy", t_data, 0, 0, 0, "dwdy");

		// Aij[2][2]->getDataOperation(t_data);
		// esio_field_write_double(h, "dwdz", t_data, 0, 0, 0, "dwdz");
	
		esio_file_close(h);
		esio_handle_finalize(h);

 
	}


	// Clear vectors and tensors created with *Assign
	vel.clear();
	Aij.clear();

	// Clean up temporary fields
	PFieldVectorDelete( grad_u );
	PFieldVectorDelete( grad_v );
	PFieldVectorDelete( grad_w );
	PFieldVectorDelete( vorticity );

	PFieldTensorDelete( Sij );
	PFieldTensorDelete( Tij );
	PFieldTensorDelete( Tji );

	delete S2;
	delete T2;

	delete Q;		
	delete Q1;

	delete [] t_data;

	
	delete u;
	delete v;
	delete w;
	delete p;

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

	MPI_Barrier(MPI_COMM_WORLD);

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
