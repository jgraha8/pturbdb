#include <iostream>
#include <fstream>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include "turbdb_field.hpp"
#include "esio/esio.h"

using namespace std;

namespace pturb_fields {

TurbDBField::TurbDBField(const string &db_conf_file, const int *db_dims) {
	// Set private variables
	this->db_conf_file_ = db_conf_file;
	memcpy(this->db_dims_, db_dims, sizeof(db_dims[0]) * FIELD_NDIMS);

	// Read the DB conf file
	this->readDBConfFile();

}

TurbDBField::TurbDBField(TurbDBField &turbdb_field) {

	this->db_conf_file_ = turbdb_field.getDBConfFile();
	memcpy(this->db_dims_, turbdb_field.getDBDims(),
			sizeof(this->db_dims_[0]) * FIELD_NDIMS);

	// Read the DB conf file
	this->readDBConfFile();
}

void TurbDBField::dbFieldInit(const int *field_offset, const int *dims,
		FieldDecomp_t field_decomp, const int *periodic, int operator_order) {
	// Set private data
	memcpy(this->field_offset_, field_offset,
			sizeof(field_offset[0]) * FIELD_NDIMS);
	// Initialize the field
	this->FieldInit(dims, field_decomp, periodic, operator_order);
}

string &TurbDBField::getDBConfFile() {
	return this->db_conf_file_;
}
int *TurbDBField::getDBDims() {
	return this->db_dims_;
}
int *TurbDBField::getFieldOffset() {
	return this->field_offset_;
}
vector<double> &TurbDBField::getDBTime() {
	return this->db_time_;
}
vector<string> &TurbDBField::getDBFileNames() {
	return this->db_file_names_;
}

/*
 * Reads the database configuration file: reads the grid file name and the time
 * step and field file names
 */
void TurbDBField::readDBConfFile() {

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
	for (int i = 0; i < this->db_time_.size(); i++) {
		cout << "time, file name : " << this->db_time_.at(i) << " "
				<< this->db_file_names_.at(i) << endl;
	}

	// set the time step
	this->db_time_nsteps_ = this->db_time_.size() - 2; // Number of available time steps
	this->db_time_step_ = this->db_time_.at(1) - this->db_time_.at(0);
	this->db_time_min_ = this->db_time_.at(1); // The time at 0 at t = -dt
	this->db_time_max_ = *(this->db_time_.end() - 2);

}

void TurbDBField::pchipInit() {

	// Have to initialize the data for the calculation of the PCHIP interpolation
	this->pchip_fd_.stencil[0] = -1;
	this->pchip_fd_.stencil[1] = 0;
	this->pchip_fd_.stencil[2] = 1;

	this->pchip_fd_.ds[0] = -this->db_time_step_;
	this->pchip_fd_.ds[1] = 0;
	this->pchip_fd_.ds[2] = this->db_time_step_;

	// Compute the coefficients
	autofd_stencil_coefs_c(3, this->pchip_fd_.stencil, this->pchip_fd_.ds, 1,
			this->pchip_fd_.coefs);
}
/*
 * Computes/evaluates the 4 Hermite basis functions used in PCHIP
 * interpolation. The value tau is normalized time within the interval
 *  the interpolation is performed.
 */
void TurbDBField::pchipComputeBasis(double tau, double hermite_basis[4]) {

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
void TurbDBField::pchipComputeWeights(double hermite_basis[4],
		double pchip_weights[4]) {

	pchip_weights[0] = hermite_basis[1] * this->pchip_fd_.coefs[0];
	pchip_weights[1] = hermite_basis[0]
			+ hermite_basis[1] * this->pchip_fd_.coefs[1]
			+ hermite_basis[3] * this->pchip_fd_.coefs[0];
	pchip_weights[2] = hermite_basis[2]
			+ hermite_basis[1] * this->pchip_fd_.coefs[2]
			+ hermite_basis[3] * this->pchip_fd_.coefs[1];
	pchip_weights[3] = hermite_basis[3] * this->pchip_fd_.coefs[2];
}

/*
 * Reads the grid file and returns the x, y, and z grid coordinates for the
 * local domain. Each process will read the entire grid file and then copy
 * the relevant portion the appropriate array.
 */
void TurbDBField::readDBGridLocal(double *x, double *y, double *z) {
	// Master process reads the entire grid to a buffer
	if (this->getMpiTopology()->rank == 0) {

	}
	// The entire grid is broadcast to all processes

	// Copy the correct portions to the appropriate array

}
void TurbDBField::readDBField(double time, const char *field_name) {

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

	int offset[3] = {
		this->field_offset_[0] + offset_local[0] + offset_operation[0],
		this->field_offset_[1] + offset_local[1] + offset_operation[1],
		this->field_offset_[2] + offset_local[2] + offset_operation[2]
		};

	// We now evaluate the
	for (int i = 0; i < 4; i++) {

		esio_handle h = esio_handle_initialize(this->getMpiTopology()->comm);

		// Open the database file
		const char *db_file_name = this->db_file_names_.at(cell_index - 1 + i).c_str();
		esio_file_open(h, db_file_name, 0); // Open read-only

		esio_field_establish(h, 
				     this->db_dims_[0], offset[0], dims_operation[0], 
				     this->db_dims_[1], offset[1], dims_operation[1], 
				     this->db_dims_[2], offset[2], dims_operation[2]
				     );

		esio_field_readv_float(h, field_name, data_buffer, 0, 0, 0, 1);

		esio_file_close (h);
		esio_handle_finalize(h);

		// Apply the appropriate weights to the input data
		long j = 0;
		while (j != data_buffer_size) data_buffer[j] = pchip_weights[i] * data_buffer[j++];

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
