#include <iostream>
#include <cstdlib>
#include <complex.h>
#include <cmath>
#include "lapacke.h"
#include "clock.hpp"
#include "pfield_math.hpp"

#define ZERO 1.0e-12L

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

PFieldComplexVector_t PFieldComplexVectorNew( const PField &pfield )
{
	PFieldComplexVector_t pfield_vector(2);

	for(int i=0; i<2; i++ ) 
		pfield_vector[i] = PFieldVectorNew( pfield );

	return pfield_vector;
}

void PFieldVectorDelete( PFieldVector_t &pfield_vector )
{
	for(size_t i=0; i<pfield_vector.size(); i++ )
		delete pfield_vector[i];
	pfield_vector.clear();
}

void PFieldVectorDelete( PFieldComplexVector_t &pfield_vector )
{
	for(int i=0; i<2; i++ ) 
		PFieldVectorDelete( pfield_vector[i] );
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

PFieldComplexTensor_t PFieldComplexTensorNew( const PField &pfield )
{
	PFieldComplexTensor_t tensor(2);
	tensor[0] = PFieldTensorNew( pfield );
	tensor[1] = PFieldTensorNew( pfield );
	return tensor;
}

void PFieldTensorDelete( PFieldTensor_t &pfield_tensor )
{
	PFieldTensor_t::iterator t;
	for( t=pfield_tensor.begin(); t != pfield_tensor.end(); t++ )
		PFieldVectorDelete( *t );
	pfield_tensor.clear();
}

void PFieldTensorDelete( PFieldComplexTensor_t &tensor )
{
	PFieldTensorDelete( tensor[0] );
	PFieldTensorDelete( tensor[1] );
	tensor.clear();
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

/*
 * Computes the determinant of the provided tensor.
 */
PField &PFieldTensorDeterminant( PField &det, const PFieldTensor_t &a ) 
{

	vector<vector<const double *> > b(3);
	for(int i=0; i<3; i++ ) {
		b[i].resize(3);
		for(int j=0; j<3; j++) {
			b[i][j] = a[i][j]->getDataLocal();
		}
	}

	PFIELD_LOOP_OPERATION_TO_LOCAL(a[0][0])
		det.setDataLocal( _index, 
				  b[0][0][_index] * ( b[1][1][_index] * b[2][2][_index] - b[1][2][_index] * b[2][1][_index] ) +
				  b[0][1][_index] * ( b[1][2][_index] * b[2][0][_index] - b[1][0][_index] * b[2][2][_index] ) +
				  b[0][2][_index] * ( b[1][0][_index] * b[2][1][_index] - b[1][1][_index] * b[2][0][_index] ) );
	PFIELD_LOOP_END
	return det;
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

/*
 * Computes the eigen pair (eigenvalue, eigenvector) for the given
 * tensor. The eigenvectors are stored in row major order.
 */ 
void PFieldEigenPair( PFieldComplexVector_t &eigen_value, PFieldComplexTensor_t &eigen_vector, const PFieldTensor_t &a ) 
{

	Clock clock_compute;
	Clock clock_set_matrix;
	Clock clock_set_eigenpair;
	// Local matrix for each element in the field
	lapack_complex_double *A_elem            = (lapack_complex_double *)calloc(PFIELD_NDIMS*PFIELD_NDIMS,sizeof(lapack_complex_double));
	lapack_complex_double *eigen_value_elem  = (lapack_complex_double *)calloc(PFIELD_NDIMS,             sizeof(lapack_complex_double));
	lapack_complex_double *eigen_vector_elem = (lapack_complex_double *)calloc(PFIELD_NDIMS*PFIELD_NDIMS,sizeof(lapack_complex_double));

	if( a[0][0]->getMPITopology()->rank == 0 )
		cout << "        PFieldEigenPair times:\n";

	PFIELD_LOOP_OPERATION_TO_LOCAL(a[0][0])

        clock_set_matrix.resume();
	for( int m=0; m<PFIELD_NDIMS; m++ ) {
		for( int n=0; n<PFIELD_NDIMS; n++ ) {
			// A_elem is in column major order
			A_elem[m + n*PFIELD_NDIMS] = lapack_make_complex_double(a[m][n]->data_local[_index],0.0L);
		}
	}
	clock_set_matrix.stop();

	clock_compute.resume();
	lapack_int info = LAPACKE_zgeev(LAPACK_COL_MAJOR,'N','V', PFIELD_NDIMS, A_elem, PFIELD_NDIMS, 
					eigen_value_elem, NULL, PFIELD_NDIMS, eigen_vector_elem, PFIELD_NDIMS );
	clock_compute.stop();

	/* Check for convergence */
	if( info > 0 ) {
		printf( "PFieldEigenPair: failed to compute eigenvalues.\n" );
		int ierr=0;
		MPI_Abort( a[0][0]->getMPITopology()->comm, ierr );
	}
	
	clock_set_eigenpair.resume();
	for( int m=0; m<PFIELD_NDIMS; m++ ) {

		eigen_value[0][m]->data_local[_index] = creal(eigen_value_elem[m]);
		eigen_value[1][m]->data_local[_index] = cimag(eigen_value_elem[m]);

		for( int n=0; n<PFIELD_NDIMS; n++ ) {
						
			// Careful the eigenvectors are given as rows;
			eigen_vector[0][m][n]->data_local[_index] = creal(eigen_vector_elem[m*PFIELD_NDIMS+n]);
			eigen_vector[1][m][n]->data_local[_index] = cimag(eigen_vector_elem[m*PFIELD_NDIMS+n]); 

		}
	}
	clock_set_eigenpair.stop();

	PFIELD_LOOP_END

	if( a[0][0]->getMPITopology()->rank == 0 ) {
		cout << "            setting A: " << clock_set_matrix.time() << "(s)\n";
		cout << "            computing eigenpair: " << clock_compute.time() << "(s)\n";
		cout << "            setting eigenpair: " << clock_set_eigenpair.time() << "(s)\n";
	}
}

/*
 * Computes the eigen pair (eigenvalue, eigenvector) for the given tensor which is assumed to be the
 * velocity gradient tensor. The eigen pairs are computed for each point in the field and if the
 * eigenvalues have 1 real component and a complex conjugate pair then the eigen pair is
 * retained. Otherwise all components are set to zero. The eigenvalues are stored as:
 *
 *     \lambda_r, \lambda_{cr}, \lambda_{ci} 
 *
 * where \lambda_r is the real eigenvalue and \lambda_{cr} +/- i\lambda_{ci} form the complex
 * conjugate pair.
 *
 * The eigenvectors are stored in row major order as:
 *
 *     v_r, v_{cr}, v_{ci}
 *
 * where v_r forms the first row and so on. The vectors v_{cr} and v_{ci} are the real and imaginary
 * components, respectively, of the original eigenvector pair given as:
 *
 *     v_{cr} +/- iv_{ci}
 *
 * For more information see Chakraborty, Balachandar, and Adrian, JFM (2005)
 */ 
void PFieldEigenPairVortex( PFieldVector_t &eigen_value, PFieldTensor_t &eigen_vector, const PFieldTensor_t &a ) 
{

	PFieldComplexVector_t eigen_value_complex  = PFieldComplexVectorNew( *a[0][0] );
	PFieldComplexTensor_t eigen_vector_complex = PFieldComplexTensorNew( *a[0][0] );

	Clock clock_eigen_pair;
	Clock clock_calcs;

	if( a[0][0]->getMPITopology()->rank == 0 )
		cout << "    PFieldEigenPairVortex times: \n";

	clock_eigen_pair.start();
	// Compute the eigen pair
	PFieldEigenPair( eigen_value_complex, eigen_vector_complex, a );
	clock_eigen_pair.stop();


	clock_calcs.start();

	PFIELD_LOOP_OPERATION_TO_LOCAL(a[0][0])

	std::vector<int> real_eigval;
	std::vector<int> complex_eigval;

	for (int m=0; m<PFIELD_NDIMS; m++ ) {

		// Checking the imaginary part of the eigenvalue
		if( fabs( eigen_value_complex[1][m]->data_local[_index] ) <= ZERO ) {
			real_eigval.push_back( m );
		} else {
			complex_eigval.push_back( m );
		}
	}
	
	if( real_eigval.size() == 1 && complex_eigval.size() == 2 ) {
		// Check if they are complex conjugates

		const double e1_real = eigen_value_complex[0][complex_eigval[0]]->data_local[_index];
		const double e1_imag = eigen_value_complex[1][complex_eigval[0]]->data_local[_index];
		const double e2_real = eigen_value_complex[0][complex_eigval[1]]->data_local[_index];
		const double e2_imag = eigen_value_complex[1][complex_eigval[1]]->data_local[_index];

		if( fabs( e1_real - e2_real ) <= ZERO && 
		    fabs( e1_imag + e2_imag ) <= ZERO ) {
			// These are complex conjugates
			//std::cout << "Complex conjugates: " << e1 << ", " << e2 << std::endl;	 
			eigen_value[0]->data_local[_index] =eigen_value_complex[0][real_eigval[0]]->data_local[_index];
			eigen_value[1]->data_local[_index] =eigen_value_complex[0][complex_eigval[0]]->data_local[_index]; // Real part of the first complex eigenvalue
			eigen_value[2]->data_local[_index] =eigen_value_complex[1][complex_eigval[0]]->data_local[_index]; // Imaginary part of the firt complex eigenvalue

			for( int n=0; n<PFIELD_NDIMS; n++ ) {
				eigen_vector[0][n]->data_local[_index] = eigen_vector_complex[0][real_eigval[0]][n]->data_local[_index];
				eigen_vector[1][n]->data_local[_index] = eigen_vector_complex[0][complex_eigval[0]][n]->data_local[_index];
				eigen_vector[2][n]->data_local[_index] = eigen_vector_complex[1][complex_eigval[1]][n]->data_local[_index];
			}
			continue;
		} else {
			goto set_all_zero;
		}
		goto set_all_zero;
	}

 set_all_zero:
	for (int m=0; m<PFIELD_NDIMS; m++ ) {
		eigen_value[m]->data_local[_index] = 0.0L ;
		for (int n=0; n<PFIELD_NDIMS; n++ ) {
			eigen_vector[m][n]->data_local[_index] = 0.0L ;
		}
	}

	PFIELD_LOOP_END

	clock_calcs.stop();
	if( a[0][0]->getMPITopology()->rank == 0 ) {
		cout << "        PFieldEigenPair call: " << clock_eigen_pair.time() << "(s)\n";
		cout << "        calculations: " << clock_calcs.time() << "(s)\n";
	}

	PFieldVectorDelete( eigen_value_complex );
	PFieldTensorDelete( eigen_vector_complex );

}

}
