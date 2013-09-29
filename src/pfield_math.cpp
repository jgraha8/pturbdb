#include <iostream>
#include "pfield_math.hpp"

using namespace std;
namespace pturbdb {

PFieldVector_t PFieldVectorNew( PField &pfield )
{
	PFieldVector_t pfield_vector;

	pfield_vector.push_back( new PField( pfield, false ) );
	pfield_vector.push_back( new PField( pfield, false ) );
	pfield_vector.push_back( new PField( pfield, false ) );

	return pfield_vector;
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

	c->mul( *a[0], *b[0] ) 
		+= c->mul( *a[1], *b[1] )
		+= c->mul( *a[2], *b[2] );

	return c;
}



PFieldTensor_t PFieldTensorNew( PFieldVector_t &v1, PFieldVector_t &v2, PFieldVector_t &v3 )
{
	// Create a new tensor from a single PField object
	PFieldTensor_t pfield_tensor;

	pfield_tensor.push_back( v1 );
	pfield_tensor.push_back( v2 );
	pfield_tensor.push_back( v3 );

	return pfield_tensor;
}

/*
 * Creates a new PField tensor struct using a copy of an existing
 * field. It does not copy the field data.
 */
PFieldTensor_t PFieldTensorNew( PField &pfield ) 
{
	PFieldVector_t v1 = PFieldVectorNew( pfield );
	PFieldVector_t v2 = PFieldVectorNew( pfield );
	PFieldVector_t v3 = PFieldVectorNew( pfield );
	return PFieldTensorNew( v1, v2, v3 );
}

PFieldTensor_t PFieldTensorNew( PFieldTensor_t &tensor ) 
{
	return PFieldTensorNew( *tensor[0][0] );
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


PFieldVector_t PFieldGradient( PField &pfield )
{
	PFieldVector_t grad = PFieldVectorNew( pfield );

	grad[0]->ddx( pfield );
	grad[1]->ddy( pfield );
	grad[2]->ddz( pfield );

	return grad;
}

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

PFieldTensor_t PFieldTensorDot( PFieldTensor_t &a, PFieldTensor_t &b ) 
{
	PFieldTensor_t c = PFieldTensorNew( a );

	for( size_t i=0; i<c.size(); i++ ) {
		for (size_t j=0; j<c[i].size(); j++ ) {
			c[i][j]->mul( *a[i][0], *b[0][j] )  
				+= c[i][j]->mul( *a[i][1], *b[1][j] )
				+= c[i][j]->mul( *a[i][2], *b[2][j] );
		}
	}
	return c;
}

PField *PFieldTensorDotDot( PFieldTensor_t &a, PFieldTensor_t &b )
{
	PFieldTensor_t c = PFieldTensorDot( a, b );
	return PFieldTensorTrace( c );
}

PField *PFieldTensorTrace( PFieldTensor_t &a ) 
{
	// Create a new PField class
	PField *b = new PField( *a[0][0] );
	
	b->mul( *a[0][0], *a[0][0] ) 
		+= b->mul( *a[1][1], *a[1][1] )
		+= b->mul( *a[2][2], *a[2][2] );

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



