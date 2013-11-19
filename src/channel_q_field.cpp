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

#define FIELD_NZ 384
#define FIELD_NY 256
#define FIELD_NX 512

// #define FIELD_NZ 64
// #define FIELD_NY 64
// #define FIELD_NX 64

#define H5_OUTPUT_PATH "/datascope/tdbchannel/analysis/channel-q-field"

using namespace std;
using namespace pturbdb;

int main(int argc, char *argv[]) {
	double wtime;

	MPI_Init(&argc, &argv);

	const int db_dims[] = { DB_NZ, DB_NY, DB_NX };
	//	int db_field_offset[] = { FIELD_NZ, (DB_NY-FIELD_NY)/2, FIELD_NX };
	const int db_field_offset[] = { 0, 0, 0 };
	const int field_dims[] = { FIELD_NZ, FIELD_NY, FIELD_NX };
	const int periodic[] = { 0, 0, 0 };

	PTurbDBField *u = new PTurbDBField("turbdb.conf", db_dims);

	//u->PFieldInit(db_field_offset, field_dims, FIELD_DECOMP_SLAB, periodic, 6);
	u->PFieldInit(db_field_offset, field_dims, FIELD_DECOMP_PENCIL, periodic, 6);
	
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
	const double start_time = u->getDBTimeMax()*7.1L/8.0L;
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
	
	// PField *S2 = new PField( *u, false );
	// PField *T2 = new PField( *u, false );
	PField *Q  = new PField( *u, false );

	// Temporary data buffer
	double *t_data = new double[u->getSizeOperation()];

	Clock clock;
	Clock clock_data_read;
	Clock clock_calcs;
	Clock clock_data_write;

	double time = start_time;

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

	clock_data_read.stop();

	clock_calcs.start();

	// Filter the fields read from file
	if( mpi_topology->rank == 0 ) {
		cout << "Filtering u, v, w ... ";
	}
	clock.start();
	u->filter( 6 ); 
	v->filter( 6 ); 
	w->filter( 6 ); 
	clock.stop();
	if( mpi_topology->rank == 0 ) cout << clock.time() << "(s)\n";

	// for(int p=0; p<u->getMPITopology()->nproc; p++) {
	// 	if( p == u->getMPITopology()->rank) {
	// 		printf("================================================================================\n");
	// 		printf("Process %d\n", u->getMPITopology()->rank);
	// 		printf("dims = %d, %d, %d\n", u->getDims()[0], u->getDims()[1], u->getDims()[2]);
	// 		printf("dims_local = %d, %d, %d\n", u->getDimsLocal()[0], u->getDimsLocal()[1], u->getDimsLocal()[2]);
	// 		printf("dims_operation = %d, %d, %d\n", u->getDimsOperation()[0], u->getDimsOperation()[1], u->getDimsOperation()[2]);
	// 		printf("offset_local = %d, %d, %d\n", u->getOffsetLocal()[0], u->getOffsetLocal()[1], u->getOffsetLocal()[2]);
	// 		printf("offset_operation = %d, %d, %d\n", u->getOffsetOperation()[0], u->getOffsetOperation()[1], u->getOffsetOperation()[2]);
	// 		printf("rind_size = %d\n", u->getRindSize());
	// 	}
	// 	MPI_Barrier(u->getMPITopology()->comm);
	// }
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

	PFieldTensorAssign( Aij, grad_u, grad_v, grad_w );

	MPI_Barrier( mpi_topology->comm );

	if( mpi_topology->rank == 0 ) cout << "Computing Q invariant ... ";
	clock.start();
	// Q->sub( *T2, *S2 ) *= 0.5;
	PFieldTensorDotDot( *Q, Aij, Aij ) *= (double)-0.5L; // Apply the -1/2 factor
	clock.stop();
	if( mpi_topology->rank == 0 ) cout << clock.time() << "(s)\n";

	clock_calcs.stop();

	clock_data_write.start();

	if( mpi_topology->rank == 0 ) cout << "Writing output" << endl;
		
	// Output data
	esio_handle h = esio_handle_initialize(mpi_topology->comm);

	//std::stringstream *h5file = new std::stringstream;
	// *h5file << "/datascope/tdbchannel/analysis/q-wall-time-";
	// h5file->precision(8);
	// *h5file << time << ".h5";
		
 	static const std::string h5file = std::string(H5_OUTPUT_PATH) + std::string("/q-fd6-filtered-6-pencil.h5");

	// Open the database file
	esio_file_create(h, h5file.c_str(), 1);
		
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

	Q->getDataOperation(t_data);
	esio_field_write_double(h, "Q", t_data, 0, 0, 0, "Q");

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

	// Clear vectors and tensors created with *Assign
	vel.clear();
	Aij.clear();

	// Clean up temporary fields
	PFieldVectorDelete( grad_u );
	PFieldVectorDelete( grad_v );
	PFieldVectorDelete( grad_w );
	PFieldVectorDelete( vorticity );

	delete [] t_data;

	delete u;
	delete v;
	delete w;
	delete Q;		

	MPI_Barrier(MPI_COMM_WORLD);

	// Terminate MPI.
	MPI_Finalize();

	return 0;
}
