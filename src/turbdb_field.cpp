#include <iostream>
#include <fstream>
#include <cstdio>
#include <cstdlib>
#include "turbdb_field.hpp"

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
	this->db_time_step_ = this->db_time_.at(1) - this->db_time_.at(0);

}

void TurbDBField::pchipInit() {

	// Have to initialize the data for the calculation of the PCHIP interpolation
	this->pchip_.stencil[0] = -1;
	this->pchip_.stencil[1] = 0;
	this->pchip_.stencil[2] = 1;

	this->pchip_.ds[0] = -this->db_time_step_;
	this->pchip_.ds[1] = 0;
	this->pchip_.ds[2] = this->db_time_step_;

	// Compute the coefficients
	autofd_stencil_coefs_c(3, this->pchip_.stencil, this->pchip_.ds, 1,
			this->pchip_.coef);
}

void TurbDBField::pchipComputeBasis( double tau, double hermite_coef[4] ) {

	double c1 = 1.0-tau;
	double c2 = pow(c1,2);
	double c3 = pow(tau,2);

	hermite_coef[0] = (1.0+2.0*tau)*c2;
	hermite_coef[1] = tau * c2;
	hermite_coef[2] = c3*(3.0-2.0*tau);
	hermite_coef[3] = -c3*c1;
}

/*
 * Reads the grid file and returns the x, y, and z grid coordinates for the
 * local domain. Each process will read the entire grid file and then copy
 * the relevant portion the appropriate array.
 */
void TurbDBField::readDBGridLocal(double *x, double *y, double *z) {
	// Master process reads the entire grid to a buffer
	if( this->mpi_topology_->rank == 0 ) {

	}
	// The entire grid is broadcast to all processes

	// Copy the correct portions to the appropriate array

}
void TurbDBField::readDBField(double time, const char *field_name) {
	// Create buffer the size of the operations domain


}


}
