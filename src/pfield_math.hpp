#include <vector>
#include "pfield.hpp"

namespace pturbdb {

#ifndef PFIELD_MATH_HPP_
#define PFIELD_MATH_HPP_

typedef enum {
	PFIELD_TENSOR_FULL,
	PFIELD_TENSOR_SYMMETRIC,
	PFIELD_TENSOR_ANTISYMMETRIC
} PFieldTensorType_t;

typedef std::vector<PField *> PFieldVector_t;
typedef std::vector<PFieldVector_t> PFieldTensor_t;

#endif

PFieldVector_t PFieldVectorNew( PField &pfield );
void PFieldVectorDelete( PFieldVector_t &pfield_vector );

PFieldTensor_t PFieldTensorNew( PField &pfield );
PFieldTensor_t PFieldTensorNew( PFieldVector_t &v1, PFieldVector_t &v2, PFieldVector_t &v3 );
void PFieldTensorDelete( PFieldTensor_t &pfield_tensor );

PField PFieldVectorDot( PFieldVector_t &a, PFieldVector_t &b );
PFieldVector_t PFieldGradient( PField &pfield );

PFieldTensor_t PFieldTensorSymmetric( PFieldTensor_t &tensor );
PFieldTensor_t PFieldTensorAntiSymmetric( PFieldTensor_t &tensor );
PFieldTensor_t PFieldTensorDot( PFieldTensor_t &a, PFieldTensor_t &b );
PField PFieldTensorDotDot( PFieldTensor_t &a, PFieldTensor_t &b );
PField PFieldTensorTrace( PFieldTensor_t &a );

PFieldTensor_t &PFieldTensorAdd( PFieldTensor_t &a, double b );
PFieldTensor_t &PFieldTensorSub( PFieldTensor_t &a, double b );
PFieldTensor_t &PFieldTensorMul( PFieldTensor_t &a, double b );
PFieldTensor_t &PFieldTensorDiv( PFieldTensor_t &a, double b );

}


