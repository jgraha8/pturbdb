/*
 * turb_db.hpp
 *
 *  Created on: Sep 17, 2013
 *      Author: jgraham
 */
#include <string>
#include <vector>
#include "pfield.hpp"
#include "autofd.h"

#define PCHIP_FD_STENCIL_SIZE 3

#ifndef PTURBDB_FIELD_HPP_
#define PTURBDB_FIELD_HPP_

using namespace std;

namespace pturbdb {

class PTurbDBField: public PField {// Be sure to call the empty constructor

private:

	string db_conf_file_;
	int    db_dims_[PFIELD_NDIMS];
	int    field_offset_[PFIELD_NDIMS];

	vector<double> db_time_;
	vector<string> db_file_names_;
	string         db_grid_file_name_;

	int    db_time_nsteps_;
	double db_time_step_; // Assuming constant time step
	double db_time_min_;
	double db_time_max_;

	typedef struct {
		int stencil[PCHIP_FD_STENCIL_SIZE];
		autofd_real_t ds[PCHIP_FD_STENCIL_SIZE];
		autofd_real_t coefs[PCHIP_FD_STENCIL_SIZE];  // Coefficients
	} pchip_fd_t;

	pchip_fd_t *pchip_fd_;


public:
	// Constructors
	PTurbDBField() {};
	PTurbDBField( const string &db_conf_file, const int *db_dims );
	PTurbDBField( PTurbDBField &g); // Copy constructor
	PTurbDBField( PTurbDBField &g, bool copy_field_data ); // Copy constructor

	// Deconstructors
	~PTurbDBField();

	// Procedure to set the field_offset and initialize the PField class
	void PFieldInit( const int *field_offset, const int *dims, FieldDecomp_t field_decomp, 
			  const int *periodic, int operator_order );

	// Getters and setters
	string  getDBConfFile();
	int    *getDBDims();
	int    *getFieldOffset();

	vector<double> getDBTime();
	vector<string> getDBFileNames();
	string         getDBGridFileName();

	int             getDBTimeNsteps();
	double          getDBTimeStep();
        double          getDBTimeMin();
        double          getDBTimeMax();
	pchip_fd_t     *getPCHIPFD();

        void readDBGridLocal(const char *field_names[3], double *x, double *y, double *z);
	void readDBField(double time, const char *field_name);

	void copyDataOperation( double *data );

private:

	void PTurbDBFieldCopy( PTurbDBField &g, bool copy_field_data );
	void readDBConfFile();
	void pchipInit();
	void pchipComputeBasis( double tau, double hermite_basis[4] );
	void pchipComputeWeights( double hermite_basis[4], double pchip_weights[4] );

        void syncDBGridLocal( int dim, const double *grid_operation, double *grid );
        void setDBPeriodicGridLocal( int dim, const double *grid_operation, double *grid );
	
};
}

#endif /* PTURBDB_FIELD_HPP_ */
