#include <iostream>
#include "pfield_math.hpp"

using namespace std;
namespace pturbdb {

////////////////////////////////////////////////////////////////////////////////
/// VECTOR PROCEDURES
////////////////////////////////////////////////////////////////////////////////

PFieldVector_t PFieldVectorNew( const PField &pfield )
{
	PFieldVector_t pfield_vector(3);

	for(int i=0; i<3; i++ ) 
		pfield_vector[i] = new PField( pfield, false );

	return pfield_vector;
}

void PFieldVectorDelete( PFieldVector_t &pfield_vector )
{
	for(size_t i=0; i<pfield_vector.size(); i++ )
		delete pfield_vector[i];
	pfield_vector.clear();
}

PFieldVector_t &PFieldVectorAssign( PFieldVector_t &vector, PField *a, PField *b, PField *c ) 
{
	vector[0] = a;
	vector[1] = b;
	vector[2] = c;
	return vector;
}

PField &PFieldVectorDot( PField &dot, const PFieldVector_t &a, const PFieldVector_t &b )
{
	PField *buffer = new PField( dot, false );

	dot.mul( *a[0], *b[0] );
	dot += buffer->mul( *a[1], *b[1] );
	dot += buffer->mul( *a[2], *b[2] );

	delete buffer;

	return dot;
}

PField &PFieldVectorMag( PField &mag, const PFieldVector_t &a )
{
	PFieldVectorDot( mag, a, a );
	mag.sqrt( mag );
	return mag;
}

////////////////////////////////////////////////////////////////////////////////
/// TENSOR PROCEDURES
////////////////////////////////////////////////////////////////////////////////

/*
 * Creates a new PField tensor struct using a copy of an existing
 * field. It does not copy the field data.
 */
PFieldTensor_t PFieldTensorNew( const PField &pfield ) 
{
	PFieldTensor_t tensor(3);
	PFieldVector_t v1 = PFieldVectorNew( pfield );
	PFieldVector_t v2 = PFieldVectorNew( pfield );
	PFieldVector_t v3 = PFieldVectorNew( pfield );

	tensor[0] = v1;
 	tensor[1] = v2;
 	tensor[2] = v3;

	return tensor;
}

PFieldTensor_t PFieldTensorNew( const PFieldTensor_t &tensor ) 
{
	return PFieldTensorNew( *tensor[0][0] );
}

void PFieldTensorDelete( PFieldTensor_t &pfield_tensor )
{
	PFieldTensor_t::iterator t;
	for( t=pfield_tensor.begin(); t != pfield_tensor.end(); t++ )
		PFieldVectorDelete( *t );
	pfield_tensor.clear();
}

PFieldTensor_t &PFieldTensorAssign( PFieldTensor_t &tensor, PFieldVector_t &v1, PFieldVector_t &v2, PFieldVector_t &v3 )
{
	tensor[0] = v1;
	tensor[1] = v2;
	tensor[2] = v3;
	return tensor;
}

PFieldVector_t &PFieldVectorCurl( PFieldVector_t &curl, PFieldVector_t &a )
{
	// Need buffer PField
	PField *buffer = new PField( *a[0], false );

	curl[0]->ddy( *a[2] ) -= buffer->ddz( *a[1] );
	curl[1]->ddz( *a[0] ) -= buffer->ddx( *a[2] );
	curl[2]->ddx( *a[1] ) -= buffer->ddy( *a[0] );
	
	delete buffer;
	return curl;
}


PFieldVector_t &PFieldVectorGradient( PFieldVector_t &grad, PField &pfield )
{
	grad[0]->ddx( pfield );
	grad[1]->ddy( pfield );
	grad[2]->ddz( pfield );
	return grad;
}


/*
 * Computes the symmetric component of the provided tensor
 */
PFieldTensor_t &PFieldTensorSymmetric( PFieldTensor_t &tensor_symmetric, const PFieldTensor_t &tensor ) {

	*tensor_symmetric[0][0] = *tensor[0][0];                            // Copy field 
	tensor_symmetric[0][1]->add( *tensor[0][1], *tensor[1][0] ) *= 0.5; // Compute field
	tensor_symmetric[0][2]->add( *tensor[0][2], *tensor[2][0] ) *= 0.5; // Compute field

	// Delete and copy the pointer
	*tensor_symmetric[1][0] = *tensor_symmetric[0][1]; // Copy pointer
	*tensor_symmetric[1][1] = *tensor[1][1];                                        // Copy field
	tensor_symmetric[1][2]->add( *tensor[1][2], *tensor[2][1] ) *= 0.5;             // Compute field

	*tensor_symmetric[2][0] = *tensor_symmetric[0][2]; // Copy pointer
	*tensor_symmetric[2][1] = *tensor_symmetric[1][2]; // Copy pointer
	*tensor_symmetric[2][2] = *tensor[2][2];                                        // Copy field

	return tensor_symmetric;
}

/*
 * Computes the anti-symmetric component of the provided tensor
 */
PFieldTensor_t &PFieldTensorAntiSymmetric( PFieldTensor_t &tensor_anti, const PFieldTensor_t &tensor ) {

	// Create a new base tensor
	*tensor_anti[0][0] = 0.0;                                      // Set field 
	tensor_anti[0][1]->sub( *tensor[0][1], *tensor[1][0] ) *= 0.5; // Compute field
	tensor_anti[0][2]->sub( *tensor[0][2], *tensor[2][0] ) *= 0.5; // Compute field

	// Delete and copy the pointer
	tensor_anti[1][0]->sub( *tensor[1][0], *tensor[0][1] ) *= 0.5; // Compute field
	*tensor_anti[1][1] = 0.0;                                      // Set field
	tensor_anti[1][2]->sub( *tensor[1][2], *tensor[2][1] ) *= 0.5; // Compute field

	tensor_anti[2][0]->sub( *tensor[2][0], *tensor[0][2] ) *= 0.5; // Compute field
	tensor_anti[2][1]->sub( *tensor[2][1], *tensor[1][2] ) *= 0.5; // Compute field
	*tensor_anti[2][2] = 0.0;                                      // Set field

	return tensor_anti;

}

PFieldTensor_t &PFieldTensorTranspose( PFieldTensor_t &tensor_trans, const PFieldTensor_t &tensor ) {

	for( size_t i=0; i<tensor.size(); i++ ) {
		for (size_t j=0; j<tensor[i].size(); j++ ) {
			*tensor_trans[j][i] = *tensor[i][j];
		}
	}
	return tensor_trans;
}

/*
 * Computes the inner/dot product of two tensors.
 *
 *   c = a \dot b = ab (matrix notation)
 *   c_{ij} = a_{ik}b_{kj} = a_{ik}b_{lj}\delta_{kl}
 *
 */
PFieldTensor_t &PFieldTensorDot( PFieldTensor_t &dot, const PFieldTensor_t &a, const PFieldTensor_t &b ) 
{
	PField *buffer = new PField( *a[0][0], false );

	for( size_t i=0; i<dot.size(); i++ ) {
		for (size_t j=0; j<dot[i].size(); j++ ) {
			dot[i][j]->mul( *a[i][0], *b[0][j] );
			*dot[i][j] += buffer->mul( *a[i][1], *b[1][j] );
			*dot[i][j] += buffer->mul( *a[i][2], *b[2][j] );
		}
	}
	
	delete buffer;
	return dot;
}

/*
 * Computes the inner double dot product of two tensors defined as:
 *     c = trace(ab) => c = a_{ij}b_{ji}
 */
PField &PFieldTensorDotDot( PField &dotdot, const PFieldTensor_t &a, const PFieldTensor_t &b )
{
	PFieldTensor_t c = PFieldTensorNew( dotdot );
	PFieldTensorDot( c, a, b );
	PFieldTensorTrace( dotdot, c );
	PFieldTensorDelete( c );
	return dotdot;
}

/*
 * Computes the trace of the provided tensor. Returns a new PField
 * pointer.
 */
PField &PFieldTensorTrace( PField &trace, const PFieldTensor_t &a ) 
{
	trace = *a[0][0];
	trace += *a[1][1];
	trace += *a[2][2];
	return trace;
}

PFieldTensor_t &PFieldTensorAdd( PFieldTensor_t &a, double b )
{
	for( size_t i=0; i<a.size(); i++ ) {
		for (size_t j=0; j<a[i].size(); j++ ) {
			*a[i][j] += b;
		}
	}
	return a;
}

PFieldTensor_t &PFieldTensorSub( PFieldTensor_t &a, double b )
{
	for( size_t i=0; i<a.size(); i++ ) {
		for (size_t j=0; j<a[i].size(); j++ ) {
			*a[i][j] -= b;
		}
	}
	return a;
}

PFieldTensor_t &PFieldTensorMul( PFieldTensor_t &a, double b )
{
	for( size_t i=0; i<a.size(); i++ ) {
		for (size_t j=0; j<a[i].size(); j++ ) {
			*a[i][j] *= b;
		}
	}
	return a;
}

PFieldTensor_t &PFieldTensorDiv( PFieldTensor_t &a, double b )
{
	for( size_t i=0; i<a.size(); i++ ) {
		for (size_t j=0; j<a[i].size(); j++ ) {
			*a[i][j] /= b;
		}
	}
	return a;
}


}



