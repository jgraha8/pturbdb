#include <iostream>
#include <fstream>
#include <cstdio>
#include <cstdlib>
#include "turb_db.hpp"

using namespace std;

namespace pturb_fields {

TurbDB::TurbDB(const string &db_conf_file, const int *db_dims) {
	// Set private variables
	this->db_conf_file_ = db_conf_file;
	memcpy(this->db_dims_, db_dims, sizeof(db_dims[0]) * FIELD_NDIMS);

	// Read the DB conf file
	this->readDBConfFile();

}

void TurbDB::dbFieldInit(const int *field_offset, const int *dims,
		FieldDecomp_t field_decomp, const int *periodic, int operator_order) {
	// Set private data
	memcpy(this->field_offset_, field_offset,
			sizeof(field_offset[0]) * FIELD_NDIMS);
	// Initialize the field
	this->FieldInit(dims, field_decomp, periodic, operator_order);
}

void TurbDB::readDBConfFile() {

	char db_field[80];
	double time;

	FILE *db_conf_file = fopen(this->db_conf_file_.c_str(), "r");

	if( db_conf_file == NULL ) {
		cerr << "unable to open database config file" << endl;
		exit(EXIT_FAILURE);
	}

	while( fscanf( db_conf_file, "%lf %s", &time, db_field) != EOF ) {
		this->db_time_.push_back(time);
		this->db_file_names_.push_back(db_field);
	}
	fclose(db_conf_file);

	// Print the time and database file names
	for(int i=0; i<this->db_time_.size(); i++ ) {
		cout << "time, file name : " << this->db_time_.at(i) << " " << this->db_file_names_.at(i) << endl;
	}
}
}
