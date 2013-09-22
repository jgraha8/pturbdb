#include <iostream>
#include <fstream>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include "pturbdb_field.hpp"
#include "esio/esio.h"

using namespace std;

namespace pturbdb {

PTurbDBField::PTurbDBField(const string &db_conf_file, const int *db_dims) {

	// Set private variables
	this->db_conf_file_ = db_conf_file;
	memcpy(this->db_dims_, db_dims, sizeof(db_dims[0]) * PFIELD_NDIMS);

	// Read the DB conf file
	this->readDBConfFile();

	// Initialize the PCHIP finite difference struct
	this->pchipInit();

}

/*
 * Copy Constructor
 */
PTurbDBField::PTurbDBField( PTurbDBField &g )
{
	this->PTurbDBFieldCopy( g, true );
}

PTurbDBField::PTurbDBField( PTurbDBField &g, bool copy_field_data )
{
	this->PTurbDBFieldCopy( g, copy_field_data );
}

PTurbDBField::~PTurbDBField() {
	delete this->pchip_fd_;
}

/*
 * Worker procedure for the copy constructors
 */ 
void PTurbDBField::PTurbDBFieldCopy( PTurbDBField &g, bool copy_field_data ) 
{

	// Copy/set data members in the PTurbDBField class
	this->db_conf_file_ = g.getDBConfFile();

	memcpy(this->db_dims_,      g.getDBDims(),      sizeof(this->db_dims_[0])*PFIELD_NDIMS);
	memcpy(this->field_offset_, g.getFieldOffset(), sizeof(this->field_offset_[0])*PFIELD_NDIMS);

	this->db_time_           = g.getDBTime();
	this->db_file_names_     = g.getDBFileNames();
	this->db_grid_file_name_ = g.getDBGridFileName();

	this->db_time_nsteps_ = g.getDBTimeNsteps();
	this->db_time_step_   = g.getDBTimeStep();
	this->db_time_min_    = g.getDBTimeMin();
	this->db_time_max_    = g.getDBTimeMax();

	this->pchip_fd_ = g.getPCHIPFD();

	// Copy the field class members
	this->PFieldCopy( g, copy_field_data );

}

void PTurbDBField::dbFieldInit(const int *field_offset, const int *dims,
			      FieldDecomp_t field_decomp, const int *periodic, int operator_order) {
	// Set private data
	memcpy(this->field_offset_, field_offset,
	       sizeof(field_offset[0]) * PFIELD_NDIMS);
	// Initialize the field
	this->PFieldInit(dims, field_decomp, periodic, operator_order);
}

/*
 * Data member getters
 */
string                   PTurbDBField::getDBConfFile()     { return this->db_conf_file_;      }
int                     *PTurbDBField::getDBDims()         { return this->db_dims_;           }
int                     *PTurbDBField::getFieldOffset()    { return this->field_offset_;      }
vector<double>           PTurbDBField::getDBTime()         { return this->db_time_;           }
vector<string>           PTurbDBField::getDBFileNames()    { return this->db_file_names_;     }
string                   PTurbDBField::getDBGridFileName() { return this->db_grid_file_name_; }
int                      PTurbDBField::getDBTimeNsteps()   { return this->db_time_nsteps_;    }
double                   PTurbDBField::getDBTimeStep()     { return this->db_time_step_;      }
double                   PTurbDBField::getDBTimeMin()      { return this->db_time_min_;       }
double                   PTurbDBField::getDBTimeMax()      { return this->db_time_max_;       }
PTurbDBField::pchip_fd_t *PTurbDBField::getPCHIPFD()        { return this->pchip_fd_;          }

/*
 * Reads the database configuration file: reads the grid file name and the time
 * step and field file names
 */
void PTurbDBField::readDBConfFile() {

	char db_field[80];
	double time;

	FILE *db_conf_file = fopen(this->db_conf_file_.c_str(), "r");

	if (db_conf_file == NULL) {
		cerr << "unable to open database config file" << endl;
		exit(EXIT_FAILURE);
	}

	// The first line has the grid file
	if (fscanf(db_conf_file, "%s", db_field) != EOF) {
		this->db_grid_file_name_ = db_field;
	}
	while (fscanf(db_conf_file, "%lf %s", &time, db_field) != EOF) {
		this->db_time_.push_back(time);
		this->db_file_names_.push_back(db_field);
	}
	fclose(db_conf_file);

	// Print the time and database file names
	// for (int i = 0; i < this->db_time_.size(); i++) {
	// 	cout << "time, file name : " << this->db_time_.at(i) << " "
	// 			<< this->db_file_names_.at(i) << endl;
	// }

	// set the time step
	this->db_time_nsteps_ = this->db_time_.size() - 2; // Number of available time steps
	this->db_time_step_   = this->db_time_.at(1) - this->db_time_.at(0);
	this->db_time_min_    = this->db_time_.at(1); // The time at 0 at t = -dt
	this->db_time_max_    = *(this->db_time_.end() - 2);

}

void PTurbDBField::pchipInit() {

	// Create new memory block
	this->pchip_fd_ = new pchip_fd_t;

	// Have to initialize the data for the calculation of the PCHIP interpolation
	this->pchip_fd_->stencil[0] = -1;
	this->pchip_fd_->stencil[1] = 0;
	this->pchip_fd_->stencil[2] = 1;

	this->pchip_fd_->ds[0] = -this->db_time_step_;
	this->pchip_fd_->ds[1] = 0;
	this->pchip_fd_->ds[2] = this->db_time_step_;

	// Compute the coefficients
	autofd_stencil_coefs_c(3, this->pchip_fd_->stencil, this->pchip_fd_->ds, 1,
			       this->pchip_fd_->coefs);
}
/*
 * Computes/evaluates the 4 Hermite basis functions used in PCHIP
 * interpolation. The value tau is normalized time within the interval
 *  the interpolation is performed.
 */
void PTurbDBField::pchipComputeBasis(double tau, double hermite_basis[4]) {

	double c1 = 1.0 - tau;
	double c2 = pow(c1, 2);
	double c3 = pow(tau, 2);

	hermite_basis[0] = (1.0 + 2.0 * tau) * c2;
	hermite_basis[1] = tau * c2;
	hermite_basis[2] = c3 * (3.0 - 2.0 * tau);
	hermite_basis[3] = -c3 * c1;
}
/*
 * Computes/evaluates the weights used in the PCHIP interpolation. The weights
 * w are applied as:
 *
 *   p(t) = \sum_i=0^3 w_i p_(n-1 + i)
 *
 * and the integer n determines the interval the interpolation is performed over.
 * The discrete values p_j is the known field
 */
void PTurbDBField::pchipComputeWeights(double hermite_basis[4],
				      double pchip_weights[4]) {

	if( this->pchip_fd_ == NULL ) {
		cout << "PTurbDBField: PCHIP finite difference struct not allocated" << endl;
		int ierr=0;
		MPI_Abort(this->getMpiTopology()->comm, ierr);
	}
	
	pchip_weights[0] = hermite_basis[1] * this->pchip_fd_->coefs[0];
	pchip_weights[1] = hermite_basis[0]
		+ hermite_basis[1] * this->pchip_fd_->coefs[1]
		+ hermite_basis[3] * this->pchip_fd_->coefs[0];
	pchip_weights[2] = hermite_basis[2]
		+ hermite_basis[1] * this->pchip_fd_->coefs[2]
		+ hermite_basis[3] * this->pchip_fd_->coefs[1];
	pchip_weights[3] = hermite_basis[3] * this->pchip_fd_->coefs[2];
}

/*
 * Reads the grid file and returns the x, y, and z grid coordinates for the
 * local domain. Each process will read the entire grid file and then copy
 * the relevant portion the appropriate array.
 */
void PTurbDBField::readDBGridLocal(const char *field_names[3], double *x,
				  double *y, double *z) {

	// Get the offset for the local domain
	const int *offset_local = this->getOffsetLocal();

	const int *offset_operation = this->getOffsetOperation();
	const int *dims_operation = this->getDimsOperation();

	int offset[3] = { this->field_offset_[0] + offset_local[0]
			  + offset_operation[0], this->field_offset_[1] + offset_local[1]
			  + offset_operation[1], this->field_offset_[2] + offset_local[2]
			  + offset_operation[2] };

	// Need buffers to store the grid for the operation domain
	double *x_operation = new double[dims_operation[0]];
	double *y_operation = new double[dims_operation[1]];
	double *z_operation = new double[dims_operation[2]];

	esio_handle h = esio_handle_initialize(this->getMpiTopology()->comm);

	// Open the database file
	const char *db_grid_file_name = this->db_grid_file_name_.c_str();

	esio_file_open(h, db_grid_file_name, 0); // Open read-only

	// Read x
	esio_line_establish(h, this->db_dims_[0], offset[0], dims_operation[0]);
	esio_line_read_double(h, field_names[0], x_operation, 0);
	// Read y
	esio_line_establish(h, this->db_dims_[1], offset[1], dims_operation[1]);
	esio_line_read_double(h, field_names[1], y_operation, 0);
	// Read z
	esio_line_establish(h, this->db_dims_[2], offset[2], dims_operation[2]);
	esio_line_read_double(h, field_names[2], z_operation, 0);

	esio_file_close(h);
	esio_handle_finalize(h);

	// First set the interior values of grid from grid_operation
	for (long i = 0; i < dims_operation[0]; i++)
		x[offset_operation[0] + i] = x_operation[i];
	for (long i = 0; i < dims_operation[1]; i++)
		y[offset_operation[1] + i] = y_operation[i];
	for (long i = 0; i < dims_operation[2]; i++)
		z[offset_operation[2] + i] = z_operation[i];

	if (this->getFieldDecomp() == FIELD_DECOMP_SLAB
	    || this->getFieldDecomp() == FIELD_DECOMP_PENCIL) {
		this->syncDBGridLocal(0, x_operation, x);
	}

	if (this->getFieldDecomp() == FIELD_DECOMP_PENCIL) {
		this->syncDBGridLocal(1, y_operation, y);
	}

	// Next piece if we have periodic domain and we are a boundary
	// process, we have to manually set the rind grid values
	this->setDBPeriodicGridLocal( 0, x_operation, x);
	this->setDBPeriodicGridLocal( 1, y_operation, y);
	this->setDBPeriodicGridLocal( 2, z_operation, z);

	// Set the local grid for the field; only sets the pointers
	this->setGridLocal( x, y, z );
	// Now initialize the finite difference class
	this->finiteDiffInit();

}

/*
 * Synchronizes and updates the grid for the local domain. The grid
 * must be read to the operation domain and the rind points are then
 * set by 1) synchronizing the grid and 2) for processes on the
 * periphery of the domain when using periodic boundary conditions,
 * they are set manually assuming a uniform grid distribution.
 */
void PTurbDBField::syncDBGridLocal(int dim, const double *grid_operation,
				  double *grid) {

	// Copy the correct portions to the appropriate array

	// Must synchronize the grid data
	MPI_Request prev_request, next_request;
	MPI_Status status;

	MpiTopology_t *mpi_topology = this->getMpiTopology();

	const long rind_size = this->getRindSize();
	double *prev_buffer = new double[rind_size];
	double *next_buffer = new double[rind_size];

	int *dims_operation = this->getDimsOperation();
	int *dims_local = this->getDimsLocal();

	/*
	 * First handshakes
	 */
	// Recv next
	MPI_Irecv(next_buffer, rind_size, MPI_DOUBLE,
		  mpi_topology->neighbor_next[dim], 1, mpi_topology->comm,
		  &next_request);

	// Send prev
	// Packs the buffer at the lower indicies
	if( this->hasRind(dim,-1) ) {
	        for (long i = 0; i < rind_size; i++)
		        prev_buffer[i] = grid_operation[i];
	}

	MPI_Isend(prev_buffer, rind_size, MPI_DOUBLE,
		  mpi_topology->neighbor_prev[dim], 1, mpi_topology->comm,
		  &prev_request);

	MPI_Wait(&prev_request, &status);
	MPI_Wait(&next_request, &status);

	// Only unpack the data if it has a rind
	if( this->hasRind( dim, 1 ) ) {
	        for (long i = 0; i < rind_size; i++)
		        grid[dims_local[dim] - rind_size + i] = next_buffer[i];
	}

	/*
	 * Second handshakes
	 */
	// Recv prev
	MPI_Irecv(prev_buffer, rind_size, MPI_DOUBLE,
		  mpi_topology->neighbor_prev[dim], 2, mpi_topology->comm,
		  &prev_request);

	// Send next
	if( this->hasRind( dim, 1 ) ) {
	        for (long i = 0; i < rind_size; i++)
		        next_buffer[i] = grid_operation[dims_operation[dim] - rind_size + i];
	}

	MPI_Isend(next_buffer, rind_size, MPI_DOUBLE,
		  mpi_topology->neighbor_next[dim], 2, mpi_topology->comm,
		  &next_request);

	MPI_Wait(&next_request, &status);
	MPI_Wait(&prev_request, &status);

	// Only unpack the data if it has a rind
	if( this->hasRind( dim, -1 ) ) {
	        for (long i = 0; i < rind_size; i++)
		        grid[i] = prev_buffer[i];
	}

	delete[] prev_buffer;
	delete[] next_buffer;
}

void PTurbDBField::setDBPeriodicGridLocal( int dim, const double *grid_operation, double *grid ) {

	int *dims_local = this->getDimsLocal();
	int *dims_operation = this->getDimsOperation();

	MpiTopology_t *mpi_topology = this->getMpiTopology();
	int rind_size = this->getRindSize();

	if (this->getFieldPeriodic()[dim] == 1) {
		// Assuming uniform grid spacing when we have a periodic direction
		double ds = grid_operation[1] - grid_operation[0];

		if (mpi_topology->coords[dim] == 0) {
			for (long i = 0; i < rind_size; i++)
				grid[i] = grid_operation[0] - (rind_size - i) * ds;

		}
		if (mpi_topology->coords[dim] == mpi_topology->dims[dim] - 1) {
			for (long i = 0; i < rind_size; i++)
				grid[dims_local[dim] - rind_size + i] =
					grid_operation[dims_operation[dim] - 1] + (i + 1) * ds;

		}
	}


}

void PTurbDBField::readDBField(double time, const char *field_name) {

	// Now compute which
	if (time < this->db_time_min_ || time > this->db_time_max_) {
		cerr << "time out of bounds" << endl;
		exit(EXIT_FAILURE);
	}

	int cell_index = floor((time - this->db_time_min_) / this->db_time_step_) + 1;
	// Make sure we don't pick the last point as the cell value; this only happens when time = db_time_max_
	cell_index = fmin(cell_index, this->db_time_nsteps_ - 1);

	double tau = (time - this->db_time_.at(cell_index)) / this->db_time_step_; // 0<= tau <= 1

	// Compute the Hermite basis functions for the given normalized time.
	double hermite_basis[4], pchip_weights[4];

	this->pchipComputeBasis(tau, hermite_basis);
	this->pchipComputeWeights(hermite_basis, pchip_weights);

	// Create buffer the size of the operations domain
	const long data_buffer_size = this->getSizeOperation();
	float *data_buffer = new float[data_buffer_size];

	// Get the offset for the local and operation domains
	const int *offset_local = this->getOffsetLocal();
	const int *offset_operation = this->getOffsetOperation();
	const int *dims_operation = this->getDimsOperation();

	int offset[3] = { this->field_offset_[0] + offset_local[0]
			  + offset_operation[0], this->field_offset_[1] + offset_local[1]
			  + offset_operation[1], this->field_offset_[2] + offset_local[2]
			  + offset_operation[2] };

	// We now evaluate the
	for (int i = 0; i < 4; i++) {

		esio_handle h = esio_handle_initialize(this->getMpiTopology()->comm);

		// Open the database file
		const char *db_file_name =
			this->db_file_names_.at(cell_index - 1 + i).c_str();
		esio_file_open(h, db_file_name, 0); // Open read-only

		esio_field_establish(h, this->db_dims_[0], offset[0], dims_operation[0],
				     this->db_dims_[1], offset[1], dims_operation[1],
				     this->db_dims_[2], offset[2], dims_operation[2]);

		esio_field_read_float(h, field_name, data_buffer, 0, 0, 0);

		esio_file_close(h);
		esio_handle_finalize(h);

		// Apply the appropriate weights to the input data
		long j = 0;
		while (j != data_buffer_size) {
			data_buffer[j] = pchip_weights[i] * data_buffer[j];
			j++;
		}

		// We now modify the data_local field using the pchip_weights and the data read
		if (i == 0) {
			this->setDataOperation(data_buffer);
		} else {
			this->addDataOperation(data_buffer);
		}

	}

	// Make sure we set the synchronized_ flag to false
	this->setSynchronized(false);
}

}
