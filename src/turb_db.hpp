/*
 * turb_db.hpp
 *
 *  Created on: Sep 17, 2013
 *      Author: jgraham
 */
#include <string>
#include <vector>
#include "field.hpp"

#ifndef TURB_DB_HPP_
#define TURB_DB_HPP_

using namespace std;

namespace pturb_fields {

class TurbDB: public Field // Be sure to call the empty constructor
{

private:

	string db_conf_file_;
	int db_dims_[FIELD_NDIMS];
	int field_offset_[FIELD_NDIMS];

	vector<double> db_time_;
	vector<string> db_file_names_;

public:
	// Constructors
	TurbDB() {};
	TurbDB( const string &db_conf_file, const int *db_dims );
	// Deconstructors
	~TurbDB() {};

	void dbFieldInit( const int *field_offset, const int *dims, FieldDecomp_t field_decomp,
			const int *periodic, int operator_order );

private:

	void readDBConfFile();

};
}

#endif /* TURB_DB_HPP_ */
