/*
 * turb_db.hpp
 *
 *  Created on: Sep 17, 2013
 *      Author: jgraham
 */

#ifndef PTURBDB_FIELD_HPP_
#define PTURBDB_FIELD_HPP_

#include <string>
#include <vector>
#include "field_cache.hpp"
#include "pfield.hpp"
#include "autofd.h"

#define PCHIP_FD_STENCIL_SIZE 3

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
	
	// Data caching mechanism
	bool pchip_caching_;
	FieldCache<float> *pchip_data_cache_;

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
	string &getDBConfFile()            { return this->db_conf_file_;      }
	int    *getDBDims()                { return this->db_dims_;           }
	int    *getFieldOffset()           { return this->field_offset_;      }

	vector<double> &getDBTime()        { return this->db_time_;           }
	vector<string> &getDBFileNames()   { return this->db_file_names_;     }
	string         &getDBGridFileName(){ return this->db_grid_file_name_; }

	int            &getDBTimeNsteps()  { return this->db_time_nsteps_;    }
	double         &getDBTimeStep()    { return this->db_time_step_;      }
        double         &getDBTimeMin()     { return this->db_time_min_;       }
        double         &getDBTimeMax()     { return this->db_time_max_;       }
	pchip_fd_t     *getPCHIPFD()       { return this->pchip_fd_;          }
	bool           &getPCHIPCaching()  { return this->pchip_caching_;}

	void setPCHIPCaching(bool caching){this->pchip_caching_ = caching;}

        void readDBGridLocal(const char *field_names[3], double *x, double *y, double *z);
	void readDBField(double time, const char *field_name);

	void copyDataOperation( double *data );
	
	// Assignment operators; these are not inherited. We are using
	// the PField assignment operator to set data_local;
	PTurbDBField &operator=( float a ){ PField::operator=(a); return *this; }
	PTurbDBField &operator=( double a ){ PField::operator=(a); return *this; }
	PTurbDBField &operator=( const float *a ){ PField::operator=(a); return *this; }
	PTurbDBField &operator=( const double *a ){ PField::operator=(a); return *this; }
	PTurbDBField &operator=( const PTurbDBField &a ){ PField::operator=(a); return *this; }

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
