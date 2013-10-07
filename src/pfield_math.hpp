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

typedef std::vector<PFieldVector_t> PFieldComplexVector_t;
typedef std::vector<PFieldTensor_t> PFieldComplexTensor_t;

#endif

PFieldVector_t PFieldVectorNew( const PField &pfield );
PFieldComplexVector_t PFieldComplexVectorNew( const PField &pfield );
void PFieldVectorDelete( PFieldVector_t &pfield_vector );
void PFieldVectorDelete( PFieldComplexVector_t &pfield_vector );

PFieldVector_t &PFieldVectorAssign( PFieldVector_t &vector, PField *a, PField *b, PField *c );

PFieldTensor_t PFieldTensorNew( const PField &pfield );
PFieldTensor_t PFieldTensorNew( const PFieldTensor_t &a );
PFieldComplexTensor_t PFieldComplexTensorNew( const PField &pfield );

void PFieldTensorDelete( PFieldTensor_t &pfield_tensor );
void PFieldTensorDelete( PFieldComplexTensor_t &pfield_tensor );

PFieldTensor_t &PFieldTensorAssign( PFieldTensor_t &tensor, PFieldVector_t &v1, PFieldVector_t &v2, PFieldVector_t &v3 );

PField &PFieldVectorDot( PField &dot, const PFieldVector_t &a, const PFieldVector_t &b );
PFieldVector_t &PFieldVectorCurl( PFieldVector_t &curl, PFieldVector_t &a );
PFieldVector_t &PFieldVectorGradient( PFieldVector_t &grad, PField &pfield );

PFieldTensor_t &PFieldTensorSymmetric( PFieldTensor_t &tensor_symmetric, const PFieldTensor_t &tensor );
PFieldTensor_t &PFieldTensorAntiSymmetric( PFieldTensor_t &tensor_anti, const PFieldTensor_t &tensor );
PFieldTensor_t &PFieldTensorTranspose( PFieldTensor_t &tensor_trans, const PFieldTensor_t &tensor );
PFieldTensor_t &PFieldTensorDot( PFieldTensor_t &dot, const PFieldTensor_t &a, const PFieldTensor_t &b );
PField &PFieldTensorDotDot( PField &dotdot, const PFieldTensor_t &a, const PFieldTensor_t &b );
PField &PFieldTensorTrace( PField &trace, const PFieldTensor_t &a );
PField &PFieldTensorDeterminant( PField &det, const PFieldTensor_t &a );

PFieldTensor_t &PFieldTensorAdd( PFieldTensor_t &a, double b );
PFieldTensor_t &PFieldTensorSub( PFieldTensor_t &a, double b );
PFieldTensor_t &PFieldTensorMul( PFieldTensor_t &a, double b );
PFieldTensor_t &PFieldTensorDiv( PFieldTensor_t &a, double b );

void PFieldEigenPair( PFieldComplexVector_t &eigval, PFieldComplexTensor_t &eigvec, const PFieldTensor_t &a );
void PFieldEigenPairVortex( PFieldVector_t &eigval, PFieldTensor_t &eigvec, const PFieldTensor_t &a );

}



