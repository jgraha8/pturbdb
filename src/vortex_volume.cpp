#include <iostream>
#include <sstream>
//#include <stdlib.h>
#include <cassert>
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

#define NSTEPS 100
#define FILTER_WIDTH 4 

#define Q_CRITERION 25.0
#define VORTEX_VOLUME_BIN_WIDTH 0.01
#define VORTEX_VOLUME_NBIN 1000 // starting number of bins

#define U_TAU 0.0499
#define LENGTH_VISC 0.0010020

#define H5_OUTPUT_PATH "/datascope/tdbchannel/analysis/filtered-vortex-field"

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

	u->PFieldInit(db_field_offset, field_dims, FIELD_DECOMP_SLAB, periodic, 6);
	
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
	const double start_time = u->getDBTimeMax()*7.9L/8.0L;
	static const double dt = 4*DB_DT;

	PTurbDBField *v = new PTurbDBField( *u, false ); // Do not copy u field data
	PTurbDBField *w = new PTurbDBField( *u, false ); // Do not copy u field data

	PFieldVector_t vel    = PFieldVectorNew( *u );
	PFieldVector_t grad_u = PFieldVectorNew( *u );
	PFieldVector_t grad_v = PFieldVectorNew( *u );
	PFieldVector_t grad_w = PFieldVectorNew( *u );

	PFieldTensor_t Aij = PFieldTensorNew( *u );

	
	PField *Q  = new PField( *u, false );
	PField *R  = new PField( *u, false );

	// PDF (only works with single process)
	std::vector<size_t> hist_vortex_volume(VORTEX_VOLUME_NBIN);

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

		// Filter the fields read from file
		if( mpi_topology->rank == 0 ) {
			cout << "Filtering u, v, w ... ";
		}
		clock.start();
		u->filter( FILTER_WIDTH ); 
		v->filter( FILTER_WIDTH ); 
		w->filter( FILTER_WIDTH ); 
		clock.stop();
		if( mpi_topology->rank == 0 ) cout << clock.time() << "(s)\n";

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

		VortexMap_t vortex_points;
		VortexRegionMap_t vortex_region;

		MPI_Barrier( mpi_topology->comm );

		if( mpi_topology->rank == 0 ) cout << "Computing Q invariant ... ";
		clock.start();
		// Q->sub( *T2, *S2 ) *= 0.5;
		PFieldTensorDotDot( *Q, Aij, Aij ) *= -(double)0.5L; // Apply the -1/2 factor
		clock.stop();
		if( mpi_topology->rank == 0 ) cout << clock.time() << "(s)\n";

		if( mpi_topology->rank == 0 ) cout << "Computing R invariant ... ";
		clock.start();
		PFieldTensorDeterminant( *R, Aij ) *= -(double)1.0L; // Apply the -1 factor
		clock.stop();
		if( mpi_topology->rank == 0 ) cout << clock.time() << "(s)\n";

		if( mpi_topology->rank == 0 ) cout << "Finding vortex points ... \n";
		clock.start();
		//VortexSearch( vortex_points, lambda, 5.0L, 1.0/sqrt(3) );
		VortexSearchQ( vortex_points, *Q, Q_CRITERION );
		MPI_Barrier( mpi_topology->comm );
		clock.stop();
		if( mpi_topology->rank == 0 ) cout << "    total: " << clock.time() << "(s)\n";

		if( mpi_topology->rank == 0 ) cout << "Assembling vortex regions ... \n";
		clock.start();
		VortexRegionSearch( vortex_region, vortex_points, *u );
		MPI_Barrier( mpi_topology->comm );
		clock.stop();
		if( mpi_topology->rank == 0 ) cout << "    total: " << clock.time() << "(s)\n";

		MPI_Barrier(mpi_topology->comm);

		// Vortex volume pdf
		for( VortexRegionMap_t::iterator vr = vortex_region.begin(); vr != vortex_region.end(); vr++ ) {
			// Normalize the volume by the viscous length scale
			//vr->second.volume_global /= pow( 0.0010020, 3 );
			size_t _n = floor( vr->second.volume_global / VORTEX_VOLUME_BIN_WIDTH );

			if( hist_vortex_volume.size() - 1 < _n ) 
				hist_vortex_volume.resize( _n + 1, 0 );
			
			hist_vortex_volume.at(_n)++;

		}
		
		clock_calcs.stop();

	}

	// Finish the vortex region volume pdf
	std::vector<double> pdf_vortex_volume( hist_vortex_volume.size() );
	std::vector<double> pdf_abscissa( hist_vortex_volume.size() );

	size_t nsamples = 0;
	for( std::vector<size_t>::iterator h=hist_vortex_volume.begin(); h != hist_vortex_volume.end(); h++ ) {
		nsamples += *h;
	}

	double pdf_integ=0.0L;
	for( size_t n=0; n<hist_vortex_volume.size(); n++ ) {
		pdf_vortex_volume.at(n) = hist_vortex_volume.at(n) / ( (double)nsamples * VORTEX_VOLUME_BIN_WIDTH ) ;
		pdf_abscissa.at(n) = ((double)n + (double)0.5L) * VORTEX_VOLUME_BIN_WIDTH;
		pdf_integ+=pdf_vortex_volume.at(n) * VORTEX_VOLUME_BIN_WIDTH;
	}
	printf("pdf_integ = %15.7e\n", pdf_integ);

	FILE *f = fopen("vortex_volume_pdf.txt","w");
	fprintf(f,"FIELD_NX %d\n", FIELD_NX);
	fprintf(f,"FIELD_NY %d\n", FIELD_NY);
	fprintf(f,"FIELD_NZ %d\n", FIELD_NZ);
	fprintf(f,"NSTEPS %d\n", NSTEPS);
	fprintf(f,"FILTER_WIDTH %d\n", FILTER_WIDTH);
	fprintf(f,"Q_CRITERION %15.7e\n", Q_CRITERION);
	fprintf(f,"VORTEX_VOLUME_BIN_WIDTH %15.7e\n", VORTEX_VOLUME_BIN_WIDTH); 
        fprintf(f,"VORTEX_VOLUME_NBIN %zd\n", pdf_vortex_volume.size() );
	fprintf(f,"U_TAU %15.7e\n", U_TAU);
	fprintf(f,"LENGTH_VISC %15.7e\n", LENGTH_VISC);
	for( size_t n=0; n<hist_vortex_volume.size(); n++ ) {
		fprintf(f,"%15.7e\t%15.7e\n", pdf_abscissa.at(n), pdf_vortex_volume.at(n) );
	}
	fclose(f);

	// Clear vectors and tensors created with *Assign
	vel.clear();
	Aij.clear();

	// Clean up temporary fields
	PFieldVectorDelete( grad_u );
	PFieldVectorDelete( grad_v );
	PFieldVectorDelete( grad_w );

	// PFieldTensorDelete( Sij );
	// PFieldTensorDelete( Tij );
	// PFieldTensorDelete( Tji );

	// delete S2;
	// delete T2;

	delete Q;		
	delete R;

	delete u;
	delete v;
	delete w;

	MPI_Barrier(MPI_COMM_WORLD);

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
