#include <iostream>
#include <sstream>
//#include <stdlib.h>
//#include <sys/types.h>
#include <sys/stat.h>
//#include <unistd.h>

#include <cassert>
#include <stdlib.h>
#include <mpi.h>
#include <cmath>
#include "esio/esio.h"
#include "clock.hpp"
#include "pdf.hpp"
#include "mpi_topology.hpp"
#include "pfield.hpp"
#include "pturbdb_field.hpp"
#include "pfield_math.hpp"
#include "vortex_tracking.hpp"

#define DB_NZ 1536
#define DB_NY 512
#define DB_NX 2048

#define DB_DT 0.0065

//#define FIELD_NZ 384
//#define FIELD_NY 256
//#define FIELD_NX 512

#define FIELD_NZ 128
#define FIELD_NY 256
#define FIELD_NX 128

#define NSTEPS 1000
#define FILTER_WIDTH 4 

//#define Q_CRITERION 5.0
//#define VORTEX_VOLUME_BIN_WIDTH 1.0e-8L
//#define VORTEX_VOLUME_Y_BIN_WIDTH 1.0e-2L
#define VORTEX_VOLUME_NBINS 1000 // starting number of bins
#define VORTEX_VOLUME_BIN_MAX 1.0
#define VORTEX_VOLUME_BIN_MIN 0.0

#define Y_MIN -1.0
#define Y_MAX 1.0
#define VORTEX_VOLUME_Y_NBINS 1000
#define VORTEX_VOLUME_Y_BIN_MAX 1.0
#define VORTEX_VOLUME_Y_BIN_MIN 0.0

#define PDF_WRITE_SKIP 10

#define KAPPA 0.4
#define U_TAU 0.0499
#define LENGTH_VISC 0.0010020

#define H5_OUTPUT_PATH "/datascope/tdbchannel/analysis/filtered-vortex-field"

#define RUN_DIRECTORY "v2"

using namespace std;
using namespace pturbdb;

//static const int N_Q_CRITERIA=8;
static const double q_criteria[] = {0.5, 1.0, 2.0, 4.0, 8.0, 16.0, 32.0, 64.0 };
//const double q_criteria[] = {4.0};
static const int N_Q_CRITERIA=sizeof(q_criteria)/sizeof(*q_criteria);

std::vector<std::vector<PDF_t> >pdf_vortex_volume(N_Q_CRITERIA);
std::vector<std::vector<PDF_2D_t> >pdf_vortex_volume_y(N_Q_CRITERIA);

void pdf_write_all();

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

	// Initialize the run directory
	struct stat st;
	if (stat(RUN_DIRECTORY, &st) == -1) {
		mkdir(RUN_DIRECTORY, 0700);
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
	static const double dt = DB_DT;

	PTurbDBField *v = new PTurbDBField( *u, false ); // Do not copy u field data
	PTurbDBField *w = new PTurbDBField( *u, false ); // Do not copy u field data

	PFieldVector_t vel    = PFieldVectorNew( *u );
	PFieldVector_t grad_u = PFieldVectorNew( *u );
	PFieldVector_t grad_v = PFieldVectorNew( *u );
	PFieldVector_t grad_w = PFieldVectorNew( *u );

	PFieldTensor_t Aij = PFieldTensorNew( *u );

	
	PField *Q  = new PField( *u, false );
	PField *R  = new PField( *u, false );

	Clock clock;
	Clock clock_data_read;
	Clock clock_calcs;
	Clock clock_data_write;

	
	const double bin_width = ( VORTEX_VOLUME_BIN_MAX - VORTEX_VOLUME_BIN_MIN ) / VORTEX_VOLUME_NBINS;
	
	for( int q=0; q<N_Q_CRITERIA; q++ ) {
		
		pdf_vortex_volume[q].resize(4);
		pdf_vortex_volume_y[q].resize(4);
		
		for( int p=0; p<4; p++ ) {

			PDFInit( pdf_vortex_volume[q][p], VORTEX_VOLUME_NBINS, bin_width, VORTEX_VOLUME_BIN_MIN );
		
			const size_t _nb[] = { VORTEX_VOLUME_Y_NBINS, VORTEX_VOLUME_NBINS };
			const double _bw[] = { ( VORTEX_VOLUME_Y_BIN_MAX - VORTEX_VOLUME_Y_BIN_MIN ) / VORTEX_VOLUME_Y_NBINS,
					       ( VORTEX_VOLUME_BIN_MAX - VORTEX_VOLUME_BIN_MIN ) / VORTEX_VOLUME_NBINS };
			const double _az[] = { VORTEX_VOLUME_Y_BIN_MIN, VORTEX_VOLUME_BIN_MIN };
			
			std::vector<size_t> _nbins(_nb, _nb + sizeof(_nb)/sizeof(*_nb));
			std::vector<double> _bin_width(_bw, _bw + sizeof(_bw)/sizeof(*_bw));
			std::vector<double> _abscissa_min(_az, _az + sizeof(_az)/sizeof(*_az));

			PDFInit( pdf_vortex_volume_y[q][p], _nbins, _bin_width, _abscissa_min );
		}
		
	}

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

		// For each Q criterion find the vortices

		if( mpi_topology->rank == 0 ) cout << "Binning data for each Q criterion ... \n";
			
		for (int q=0; q<N_Q_CRITERIA; q++ ) {
		
			if( mpi_topology->rank == 0 ) printf("    Q criterion: %15.7f (%d of %d)\n", q_criteria[q], q+1, N_Q_CRITERIA);
			if( mpi_topology->rank == 0 ) cout << "    Finding vortex points ... \n";
			clock.start();
			//VortexSearch( vortex_points, lambda, 5.0L, 1.0/sqrt(3) );
			VortexSearchQ( vortex_points, *Q, q_criteria[q] );
			MPI_Barrier( mpi_topology->comm );
			clock.stop();
			if( mpi_topology->rank == 0 ) cout << "        total: " << clock.time() << "(s)\n";
				
			if( mpi_topology->rank == 0 ) cout << "    Assembling vortex regions ... \n";
			clock.start();
			VortexRegionSearch( vortex_region, vortex_points, *Q );
			MPI_Barrier( mpi_topology->comm );
			clock.stop();
			if( mpi_topology->rank == 0 ) cout << "        total: " << clock.time() << "(s)\n";
			
			MPI_Barrier(mpi_topology->comm);
				
			if( mpi_topology->rank == 0 ) cout << "    Binning vortex data ... \n";
			clock.start();
			
			// Bin the vortex volumes
			for( VortexRegionMap_t::iterator vr = vortex_region.begin(); vr != vortex_region.end(); vr++ ) {

				const double _y = ( vr->second.barycenter[1] < 0.0 ? 
						    vr->second.barycenter[1] - Y_MIN : 
						    Y_MAX - vr->second.barycenter[1]);

				for( int p=0; p<4; p++ ) {

					const double _volume = vr->second.volume / pow( KAPPA * _y, p );

					PDFBinSample( pdf_vortex_volume[q][p], _volume );

					//printf("y barycenter = %.7e\n", vr->second.barycenter[1]);

					const double _s[] = { _y, _volume };

					std::vector<double> _sample(_s, _s + sizeof(_s)/sizeof(*_s) );
					PDFBinSample( pdf_vortex_volume_y[q][p], _sample );

				}
			}
			clock.stop();
			if( mpi_topology->rank == 0 ) cout << "        total: " << clock.time() << "(s)\n";

		}

		// Check if the pdfs should be written
		if( ( n + 1 ) % PDF_WRITE_SKIP == 0 ) 
			pdf_write_all();

		clock_calcs.stop();

	}

	// Write the final pdfs
	pdf_write_all();

	
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

void pdf_write_all()
{

  	for( int q=0; q<N_Q_CRITERIA; q++ ) {

		for( int p=0; p<4; p++ ) {
			
			PDFComputeAll( pdf_vortex_volume[q][p] );
			PDFComputeAll( pdf_vortex_volume_y[q][p] );

			char fname[80];
			sprintf(fname,"%s/vortex_volume_pdf_q_%.7f_%d.txt", RUN_DIRECTORY, q_criteria[q], p );
			FILE *f = fopen(fname,"w");
			fprintf(f,"FIELD_NX %d\n", FIELD_NX);
			fprintf(f,"FIELD_NY %d\n", FIELD_NY);
			fprintf(f,"FIELD_NZ %d\n", FIELD_NZ);
			fprintf(f,"NSTEPS %d\n", NSTEPS);
			fprintf(f,"FILTER_WIDTH %d\n", FILTER_WIDTH);
			fprintf(f,"Q_CRITERION %15.7e\n", q_criteria[q]);
			fprintf(f,"NBINS %zd\n", pdf_vortex_volume[q][p].nbins); 
			fprintf(f,"BIN_WIDTH %15.7e\n", pdf_vortex_volume[q][p].bin_width); 
			fprintf(f,"BIN_MAX %15.7e\n", (double)VORTEX_VOLUME_BIN_MAX ); 
			fprintf(f,"BIN_MIN %15.7e\n", (double)VORTEX_VOLUME_BIN_MIN ); 
			fprintf(f,"NSAMPLES %zd\n", pdf_vortex_volume[q][p].nsamples );
			fprintf(f,"U_TAU %15.7e\n", U_TAU);
			fprintf(f,"LENGTH_VISC %15.7e\n", LENGTH_VISC);
			for( size_t n=0; n<pdf_vortex_volume[q][p].nbins; n++ ) {
				fprintf(f,"%15.7e\t%15.7e\n", pdf_vortex_volume[q][p].abscissa.at(n), pdf_vortex_volume[q][p].pdf.at(n) );
			}
			fclose(f);

			sprintf(fname,"%s/vortex_volume_pdf_y_q_%.7f_%d.txt", RUN_DIRECTORY, q_criteria[q], p );
			f = fopen(fname,"w");
			fprintf(f,"VARIABLES=\"y\", \"volume\", \"PDF\"\n");
			fprintf(f,"ZONE DATAPACKING=BLOCK, i=%zd, j=%zd\n", pdf_vortex_volume_y[q][p].nbins[1], pdf_vortex_volume_y[q][p].nbins[0]);
			fprintf(f,"DT=(DOUBLE DOUBLE DOUBLE)\n");


			size_t _line=0;
			for( size_t i=0; i<pdf_vortex_volume_y[q][p].nbins[0]; i++ ) {
				for( size_t j=0; j<pdf_vortex_volume_y[q][p].nbins[1]; j++ ) {
					fprintf(f,"%15.7e", pdf_vortex_volume_y[q][p].abscissa[0][i]);
					if( ++_line % 2048 ==  0 ) fprintf(f,"\n");
				}
			}
			for( size_t i=0; i<pdf_vortex_volume_y[q][p].nbins[0]; i++ ) {
				for( size_t j=0; j<pdf_vortex_volume_y[q][p].nbins[1]; j++ ) {
					fprintf(f,"%15.7e", pdf_vortex_volume_y[q][p].abscissa[1][j]);
					if( ++_line % 2048 ==  0 ) fprintf(f,"\n");
				}
			}
			for( size_t i=0; i<pdf_vortex_volume_y[q][p].nbins[0]; i++ ) {
				for( size_t j=0; j<pdf_vortex_volume_y[q][p].nbins[1]; j++ ) {
					fprintf(f,"%15.7e", pdf_vortex_volume_y[q][p].pdf[i][j]);
					if( ++_line % 2048 ==  0 ) fprintf(f,"\n");
				}
			}
			fclose(f);
		}

	}


}
