#include <iostream>
#include <sstream>
//#include <stdlib.h>
#include <stdlib.h>
#include <mpi.h>
#include <cmath>
#include "esio/esio.h"
#include "clock.hpp"
#include "mpi_topology.hpp"
#include "pfield.hpp"
#include "pturbdb_field.hpp"
#include "pfield_math.hpp"
#include "vortex_tracking.hpp"

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
#define NSTEPS 5
#define H5_OUTPUT_PATH "/datascope/tdbchannel/analysis/test"

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
	
	// Turn on caching
	u->setPCHIPCaching(true);

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

	const int *dims_local             = u->getDimsLocal();
	const int *dims_operation         = u->getDimsOperation();
	const int *offset_local           = u->getOffsetLocal();
	const int *offset_operation       = u->getOffsetOperation();
	const MPITopology_t *mpi_topology = u->getMPITopology();

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

	PFIELD_LOOP_OPERATION( u )
	x_full[_index] = x_operation[_i];
	y_full[_index] = y_operation[_j];
	z_full[_index] = z_operation[_k];
	PFIELD_LOOP_END
	     
	// Set the starting time and the time step
	const double start_time = u->getDBTimeMax()/3.3;
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
	// PFieldTensor_t Sij = PFieldTensorNew( *u );
	// PFieldTensor_t Tij = PFieldTensorNew( *u );
	// PFieldTensor_t Tji = PFieldTensorNew( *u );
	
	PFieldVector_t lambda = PFieldVectorNew( *u );
	PFieldTensor_t eigvec = PFieldTensorNew( *u );
	std::vector<Vortex_t> vortex;

	// PField *S2 = new PField( *u, false );
	// PField *T2 = new PField( *u, false );
	PField *Q  = new PField( *u, false );
	PField *R  = new PField( *u, false );
	// PField *Q1 = new PField( *u, false );

	// Temporary data buffer
	double *t_data = new double[u->getSizeOperation()];

	Clock clock;
	Clock clock_data_read;
	Clock clock_calcs;
	Clock clock_data_write;

	for (int n=0; n<NSTEPS; n++) {

		if( mpi_topology->rank == 0 ) {
			printf("================================================================================\n");
			printf("Percent complete : %.2f\n", n*100.0/NSTEPS);
		}

		double time = start_time - dt * n;

		clock_data_read.start();

		if( mpi_topology->rank == 0 ) {
			cout << "Reading u from DB";
			if( u->getPCHIPCaching() ) {
				cout << " (cached) ... ";
			} else {
				cout << " ... ";
			}
		}
		clock.start(); 
		u->readDBField( time, "u" );
		clock.stop();
		if( mpi_topology->rank == 0 ) cout << clock.time() << "(s)\n";

		if( mpi_topology->rank == 0 ) {
			cout << "Reading v from DB";
			if( v->getPCHIPCaching() ) {
				cout << " (cached) ... ";
			} else {
				cout << " ... ";
			}
		}
		clock.start();
		v->readDBField( time, "v" );
		clock.stop();
		if( mpi_topology->rank == 0 ) cout << clock.time() << "(s)\n";

		if( mpi_topology->rank == 0 ) {
			cout << "Reading w from DB";
			if( w->getPCHIPCaching() ) {
				cout << " (cached) ... ";
			} else {
				cout << " ... ";
			}
		}
		clock.start();
		w->readDBField( time, "w" );
		clock.stop();
		if( mpi_topology->rank == 0 ) cout << clock.time() << "(s)\n";

		if( mpi_topology->rank == 0 ) {
			cout << "Reading p from DB";
			if( p->getPCHIPCaching() ) {
				cout << " (cached) ... ";
			} else {
				cout << " ... ";
			}
		}
		clock.start();
		p->readDBField( time, "p" );
		clock.stop();
		if( mpi_topology->rank == 0 ) cout << clock.time() << "(s)\n";

		clock_data_read.stop();

		clock_calcs.start();

		// Assign pointers to the velocity vector
		PFieldVectorAssign( vel, u, v, w );

		// Velocity gradients
		if( mpi_topology->rank == 0 ) cout << "Computing u gradient ... ";
		clock.start();
		PFieldVectorGradient( grad_u, *u );
		clock.stop();
		if( mpi_topology->rank == 0 ) cout << clock.time() << "(s)\n";

		if( mpi_topology->rank == 0 ) cout << "Computing v gradient ... ";
		clock.start();
		PFieldVectorGradient( grad_v, *v );		      
		clock.stop();
		if( mpi_topology->rank == 0 ) cout << clock.time() << "(s)\n";

		if( mpi_topology->rank == 0 ) cout << "Computing w gradient ... ";
		clock.start();
		PFieldVectorGradient( grad_w, *w );
		clock.stop();
		if( mpi_topology->rank == 0 ) cout << clock.time() << "(s)\n";

		if( mpi_topology->rank == 0 ) cout << "Computing vorticity ... ";
		clock.start();
		PFieldVectorCurl( vorticity, vel ) ;
		clock.stop();
		if( mpi_topology->rank == 0 ) cout << clock.time() << "(s)\n";


		PFieldTensorAssign( Aij, grad_u, grad_v, grad_w );

		// Compute the vortex eigen pair for Aij
		if( mpi_topology->rank == 0 ) cout << "Computing the vortex eigenpair of the velocity gradient tensor ... \n";
		clock.start();
		PFieldEigenPairVortex( lambda, eigvec, Aij );
		MPI_Barrier( mpi_topology->comm );
		clock.stop();
		if( mpi_topology->rank == 0 ) cout << "    total: " << clock.time() << "(s)\n";

		if( mpi_topology->rank == 0 ) cout << "Finding vortices ... \n";
		clock.start();
		VortexSearch( vortex, lambda, eigvec, 0.64, 1.0/sqrt(3) );
		MPI_Barrier( mpi_topology->comm );
		clock.stop();
		if( mpi_topology->rank == 0 ) cout << "    total: " << clock.time() << "(s)\n";

		// MPI_Barrier( mpi_topology->comm );
		// if( mpi_topology->rank == 0 ) cout << "Computing symmetric velocity gradient tensor" << endl;
		// PFieldTensorSymmetric( Sij, Aij );

		// MPI_Barrier( mpi_topology->comm );
		// if( mpi_topology->rank == 0 ) cout << "Computing anti-symmetric velocity gradient tensor" << endl;
		// PFieldTensorAntiSymmetric( Tij, Aij );
		
		// MPI_Barrier( mpi_topology->comm );
		// if( mpi_topology->rank == 0 ) cout << "Computing transpose of anti-symmetric velocity gradient tensor" << endl;
		// PFieldTensorTranspose( Tji, Tij );
		
		// MPI_Barrier( mpi_topology->comm );
		// if( mpi_topology->rank == 0 ) cout << "Computing scalar product of Sij" << endl;
		// PFieldTensorDotDot( *S2, Sij, Sij );

		// MPI_Barrier( mpi_topology->comm );
		// if( mpi_topology->rank == 0 ) cout << "Computing scalar product of Rij" << endl;
		// PFieldTensorDotDot( *T2, Tij, Tji );

		MPI_Barrier( mpi_topology->comm );

		if( mpi_topology->rank == 0 ) cout << "Computing Q invariant ... ";
		clock.start();
		// Q->sub( *T2, *S2 ) *= 0.5;
		PFieldTensorDotDot( *Q, Aij, Aij ) *= (double)-0.5L; // Apply the -1/2 factor
		clock.stop();
		if( mpi_topology->rank == 0 ) cout << clock.time() << "(s)\n";

		if( mpi_topology->rank == 0 ) cout << "Computing R invariant ... ";
		clock.start();
		PFieldTensorDeterminant( *R, Aij ) *= (double) -1.0L; // Apply the -1 factor
		clock.stop();
		if( mpi_topology->rank == 0 ) cout << clock.time() << "(s)\n";

		MPI_Barrier(mpi_topology->comm);

		clock_calcs.stop();

		clock_data_write.start();

		if( mpi_topology->rank == 0 ) cout << "Writing output" << endl;
		
		// Output data
		esio_handle h = esio_handle_initialize(mpi_topology->comm);

		//std::stringstream *h5file = new std::stringstream;
		// *h5file << "/datascope/tdbchannel/analysis/q-wall-time-";
		// h5file->precision(8);
		// *h5file << time << ".h5";
		
		std::string *h5file = new std::string(H5_OUTPUT_PATH);
		*h5file = *h5file + "/q-wall-time-";
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

		R->getDataOperation(t_data);
		esio_field_write_double(h, "R", t_data, 0, 0, 0, "R");

		// Q1->getDataOperation(t_data);
		// esio_field_write_double(h, "Q1", t_data, 0, 0, 0, "Q1");

		lambda[0]->getDataOperation(t_data);
		esio_field_write_double(h, "lambda1", t_data, 0, 0, 0, "lamba1");
		lambda[1]->getDataOperation(t_data);
		esio_field_write_double(h, "lambda2", t_data, 0, 0, 0, "lambda2");
		lambda[2]->getDataOperation(t_data);
		esio_field_write_double(h, "lambda3", t_data, 0, 0, 0, "lambda3");

		
		eigvec[0][0]->getDataOperation(t_data);
		esio_field_write_double(h, "vr1", t_data, 0, 0, 0, "vr1");
		eigvec[0][1]->getDataOperation(t_data);
		esio_field_write_double(h, "vr2", t_data, 0, 0, 0, "vr2");
		eigvec[0][2]->getDataOperation(t_data);
		esio_field_write_double(h, "vr3", t_data, 0, 0, 0, "vr3");

		eigvec[1][0]->getDataOperation(t_data);
		esio_field_write_double(h, "vcr1", t_data, 0, 0, 0, "vcr1");
		eigvec[1][1]->getDataOperation(t_data);
		esio_field_write_double(h, "vcr2", t_data, 0, 0, 0, "vcr2");
		eigvec[1][2]->getDataOperation(t_data);
		esio_field_write_double(h, "vcr3", t_data, 0, 0, 0, "vcr3");

		eigvec[2][0]->getDataOperation(t_data);
		esio_field_write_double(h, "vci1", t_data, 0, 0, 0, "vci1");
		eigvec[2][1]->getDataOperation(t_data);
		esio_field_write_double(h, "vci2", t_data, 0, 0, 0, "vci2");
		eigvec[2][2]->getDataOperation(t_data);
		esio_field_write_double(h, "vci3", t_data, 0, 0, 0, "vci3");

		std::fill_n(t_data, u->getSizeOperation(), 0.0);
		for( std::vector<Vortex_t>::iterator v=vortex.begin(); v != vortex.end(); v++ )
			t_data[v->index] = v->strength;
		esio_field_write_double(h, "vortex_strength", t_data, 0, 0, 0, "vortex_strength");

		std::fill_n(t_data, u->getSizeOperation(), -1.0);
		for( std::vector<Vortex_t>::iterator v=vortex.begin(); v != vortex.end(); v++ )
			t_data[v->index] = v->compactness;
		esio_field_write_double(h, "vortex_compactness", t_data, 0, 0, 0, "vortex_compactness");

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
		
		clock_data_write.stop();
		
		double total_time = clock_data_read.time() + clock_calcs.time() + clock_data_write.time();
		if( mpi_topology->rank == 0 ) {
			cout << endl;
			cout << "Total times:\n";
			cout << "    reads: " << clock_data_read.time() << " (" << 100*clock_data_read.time() / total_time << "%)" << endl;
			cout << "    calcs: " << clock_calcs.time() << " (" << 100*clock_calcs.time() / total_time << "%)" << endl;
			cout << "    writes: " << clock_data_write.time() << " (" << 100*clock_data_write.time() / total_time << "%)" << endl;
		}
	}


	// Clear vectors and tensors created with *Assign
	vel.clear();
	Aij.clear();

	// Clean up temporary fields
	PFieldVectorDelete( grad_u );
	PFieldVectorDelete( grad_v );
	PFieldVectorDelete( grad_w );
	PFieldVectorDelete( vorticity );

	// PFieldTensorDelete( Sij );
	// PFieldTensorDelete( Tij );
	// PFieldTensorDelete( Tji );

	PFieldVectorDelete( lambda );
	PFieldTensorDelete( eigvec );

	// delete S2;
	// delete T2;

	delete Q;		
	// delete Q1;

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
