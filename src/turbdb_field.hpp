/*
 * turb_db.hpp
 *
 *  Created on: Sep 17, 2013
 *      Author: jgraham
 */
#include <string>
#include <vector>
#include "field.hpp"
#include "autofd.h"

#define PCHIP_FD_STENCIL_SIZE 3

#ifndef TURBDB_FIELD_HPP_
#define TURBDB_FIELD_HPP_

using namespace std;

namespace pturb_fields {

class TurbDBField: public Field // Be sure to call the empty constructor
{

private:

	string db_conf_file_;
	int db_dims_[FIELD_NDIMS];
	int field_offset_[FIELD_NDIMS];

	vector<double> db_time_;
	vector<string> db_file_names_;
	string db_grid_file_name_;

	double db_time_step_; // Assuming constant time step

	typedef struct {
	    int stencil[PCHIP_FD_STENCIL_SIZE];
	    autofd_real_t ds[PCHIP_FD_STENCIL_SIZE];
	    autofd_real_t coef[PCHIP_FD_STENCIL_SIZE];  // Coefficients
	} pchip_t;

	pchip_t pchip_;


public:
	// Constructors
	TurbDBField() {};
	TurbDBField( const string &db_conf_file, const int *db_dims );
	TurbDBField( TurbDBField &turbdb_field ); // Copy constructor

	// Deconstructors
	~TurbDBField() {};

	// Getters and setters
	string &getDBConfFile();
	int *getDBDims();
	int *getFieldOffset();

	vector<double> &getDBTime();
	vector<string> &getDBFileNames();

	void dbFieldInit( const int *field_offset, const int *dims, FieldDecomp_t field_decomp,
			const int *periodic, int operator_order );
	void readDBGridLocal( double *x, double *y, double *z );
	void readDBField(double time, const char *field_name)

private:

	void readDBConfFile();
	void pchipInit();
	void pchipComputeBasis( double tau, double hermite_coef[4] )

};
}

#endif /* TURBDB_FIELD_HPP_ */
