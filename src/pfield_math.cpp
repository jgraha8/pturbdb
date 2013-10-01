#include <iostream>
#include "pfield_math.hpp"

using namespace std;
namespace pturbdb {

////////////////////////////////////////////////////////////////////////////////
/// VECTOR PROCEDURES
////////////////////////////////////////////////////////////////////////////////

PFieldVector_t PFieldVectorNew( PField &pfield )
{
	PFieldVector_t pfield_vector;

	pfield_vector.push_back( new PField( pfield, false ) );
	pfield_vector.push_back( new PField( pfield, false ) );
	pfield_vector.push_back( new PField( pfield, false ) );

	return pfield_vector;
}

PFieldVector_t PFieldVectorAssign( PField *a, PField *b, PField *c ) 
{
	PFieldVector_t vector;
	vector.push_back( a );
	vector.push_back( b );
	vector.push_back( c );
	return vector;
}

void PFieldVectorDelete( PFieldVector_t &pfield_vector )
{
	delete pfield_vector[0];
	delete pfield_vector[1];
	delete pfield_vector[2];

	pfield_vector.clear();
}

PField *PFieldVectorDot( PFieldVector_t &a, PFieldVector_t &b )
{
	PField *c = new PField( *a[0], false );
	PField *buffer = new PField( *a[0], false );

	c->mul( *a[0], *b[0] );
	*c += buffer->mul( *a[1], *b[1] );
	*c += buffer->mul( *a[2], *b[2] );

	delete buffer;
	return c;
}

PField *PFieldVectorMag( PFieldVector_t &a )
{
	PField *b = PFieldVectorDot( a, a );
	b->sqrt( *b );
	return b;
}

////////////////////////////////////////////////////////////////////////////////
/// TENSOR PROCEDURES
////////////////////////////////////////////////////////////////////////////////

/*
 * Creates a new PField tensor struct using a copy of an existing
 * field. It does not copy the field data.
 */
PFieldTensor_t PFieldTensorNew( PField &pfield ) 
{
	PFieldVector_t v1 = PFieldVectorNew( pfield );
	PFieldVector_t v2 = PFieldVectorNew( pfield );
	PFieldVector_t v3 = PFieldVectorNew( pfield );
	return PFieldTensorAssign( v1, v2, v3 );
}

PFieldTensor_t PFieldTensorNew( PFieldTensor_t &tensor ) 
{
	return PFieldTensorNew( *tensor[0][0] );
}


PFieldTensor_t PFieldTensorAssign( PFieldVector_t &v1, PFieldVector_t &v2, PFieldVector_t &v3 )
{
	// Create a new tensor from a single PField object
	PFieldTensor_t tensor;

	tensor.push_back( v1 );
	tensor.push_back( v2 );
	tensor.push_back( v3 );

	return tensor;
}

/*
 * Deletes the 
 */
void PFieldTensorDelete( PFieldTensor_t &pfield_tensor )
{
	PFieldVectorDelete( pfield_tensor[0] );
	PFieldVectorDelete( pfield_tensor[1] );
	PFieldVectorDelete( pfield_tensor[2] );

	pfield_tensor.clear();
}

PFieldVector_t PFieldCurl( PFieldVector_t &a )
{
	PFieldVector_t curl = PFieldVectorNew( *a[0] );

	// Need buffer PField
	PField *buffer = new PField( *a[0], false );

	curl[0]->ddy( *a[2] ) -= buffer->ddz( *a[1] );
	curl[1]->ddz( *a[0] ) -= buffer->ddx( *a[2] );
	curl[2]->ddx( *a[1] ) -= buffer->ddy( *a[0] );
	
	delete buffer;

	return curl;
}


PFieldVector_t PFieldGradient( PField &pfield )
{
	PFieldVector_t grad = PFieldVectorNew( pfield );

	grad[0]->ddx( pfield );
	grad[1]->ddy( pfield );
	grad[2]->ddz( pfield );

	return grad;
}

/*
 * Computes the symmetric component of the provided tensor
 */
PFieldTensor_t PFieldTensorSymmetric( PFieldTensor_t &tensor ) {

	// Create a new base tensor
	PFieldTensor_t tensor_symmetric = PFieldTensorNew( tensor );

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
PFieldTensor_t PFieldTensorAntiSymmetric( PFieldTensor_t &tensor ) {

	// Create a new base tensor
	PFieldTensor_t tensor_anti = PFieldTensorNew( tensor );

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

PFieldTensor_t PFieldTensorTranspose( PFieldTensor_t &tensor ) {

	// Create a new base tensor
	PFieldTensor_t tensor_trans = PFieldTensorNew( tensor );

	for( size_t i=0; i<tensor.size(); i++ ) {
		for (size_t j=0; j<tensor[i].size(); j++ ) {
			*tensor_trans[j][i] = *tensor[i][j];
		}
	}
	return tensor_trans;	
}
/*
 * Computes the inner/dot/tensor product of two tensors.
 */
PFieldTensor_t PFieldTensorDot( PFieldTensor_t &a, PFieldTensor_t &b ) 
{
	PFieldTensor_t c = PFieldTensorNew( a );

	PField *buffer = new PField( *a[0][0], false );

	for( size_t i=0; i<c.size(); i++ ) {
		for (size_t j=0; j<c[i].size(); j++ ) {
			c[i][j]->mul( *a[i][0], *b[0][j] );
			*c[i][j] += buffer->mul( *a[i][1], *b[1][j] );
			*c[i][j] += buffer->mul( *a[i][2], *b[2][j] );
		}
	}
	
	delete buffer;
	return c;
}

/*
 * Computes the scalar/double dot product of two tensors. Returns a
 * new PField pointer.
 */
PField *PFieldTensorDotDot( PFieldTensor_t &a, PFieldTensor_t &b )
{
	PFieldTensor_t c = PFieldTensorDot( a, b );
	return PFieldTensorTrace( c );
}

/*
 * Computes the trace of the provided tensor. Returns a new PField
 * pointer.
 */
PField *PFieldTensorTrace( PFieldTensor_t &a ) 
{
	// Create a new PField class
	PField *b = new PField( *a[0][0], false );
	PField *buffer = new PField( *a[0][0], false );

	b->mul( *a[0][0], *a[0][0] );
	*b += buffer->mul( *a[1][1], *a[1][1] );
	*b += buffer->mul( *a[2][2], *a[2][2] );

	delete buffer;
	return b;
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



