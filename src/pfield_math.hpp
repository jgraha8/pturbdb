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

PFieldVector_t &PFieldVectorAssign( PFieldVector_t &vector, PField *a, PField *b, PField *c );

PFieldTensor_t PFieldTensorNew( PField &pfield );
PFieldTensor_t PFieldTensorNew( PFieldTensor_t &a );
void PFieldTensorDelete( PFieldTensor_t &pfield_tensor );

PFieldTensor_t &PFieldTensorAssign( PFieldTensor_t &tensor, PFieldVector_t &v1, PFieldVector_t &v2, PFieldVector_t &v3 );

PField &PFieldVectorDot( PField &dot, PFieldVector_t &a, PFieldVector_t &b );
PFieldVector_t &PFieldVectorCurl( PFieldVector_t &curl, PFieldVector_t &a );
PFieldVector_t &PFieldVectorGradient( PFieldVector_t &grad, PField &pfield );

PFieldTensor_t &PFieldTensorSymmetric( PFieldTensor_t &tensor_symmetric, PFieldTensor_t &tensor );
PFieldTensor_t &PFieldTensorAntiSymmetric( PFieldTensor_t &tensor_anti, PFieldTensor_t &tensor );
PFieldTensor_t &PFieldTensorTranspose( PFieldTensor_t &tensor_trans, PFieldTensor_t &tensor );
PFieldTensor_t &PFieldTensorDot( PFieldTensor_t &dot, PFieldTensor_t &a, PFieldTensor_t &b );
PField &PFieldTensorDotDot( PField &dotdot, PFieldTensor_t &a, PFieldTensor_t &b );
PField &PFieldTensorTrace( PField &trace, PFieldTensor_t &a );

PFieldTensor_t &PFieldTensorAdd( PFieldTensor_t &a, double b );
PFieldTensor_t &PFieldTensorSub( PFieldTensor_t &a, double b );
PFieldTensor_t &PFieldTensorMul( PFieldTensor_t &a, double b );
PFieldTensor_t &PFieldTensorDiv( PFieldTensor_t &a, double b );

}



