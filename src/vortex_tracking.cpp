#include <cmath>
#include <cassert>
#include <cstdlib>
#include <cstdio>
#include <cstring>
#include <iostream>
#include <vector>
#include <deque>
#include <map>
#include <utility>
#include "lapacke.h"
#include "emap.hpp"
#include "core.hpp"
#include "pfield.hpp"
#include "pfield_math.hpp"
#include "vortex_tracking.hpp"

namespace pturbdb {

// Global variables
int ierr=0;
size_t vortex_region_tag=0;

#ifdef BOUNDS_CHECK
Vortex_t *vp_lbound, *vp_ubound;
#endif

/****************************************************************************************/
void VortexInit( Vortex_t &vortex )
/****************************************************************************************/
{
	vortex.index=-1;
	vortex.strength=0.0L;
	vortex.compactness=0.0L;
	vortex.volume = 0.0L;
}

/****************************************************************************************/
void VortexSet( Vortex_t &vortex, size_t index, double strength, double compactness, double volume ) 
/****************************************************************************************/
{
	vortex.index = index;
	vortex.strength = strength;
	vortex.compactness = compactness;
	vortex.volume = volume;
}

// void VortexCopy( Vortex_t &a, const Vortex_t &b )
// {
// 	VortexSet( a, b.index, b.strength, b.compactness );
// }

void VortexSearchQ( VortexMap_t &vortex, const PField &Q, double q_threshold ) 
{

	PFieldVector_t lambda = PFieldVectorNew( Q );

	*(lambda[0]) = (double)0.0L; // lambda_r
	*(lambda[1]) = (double)0.0L; // lambda_cr
	*(lambda[2]) = Q;    // lambda_ci; vortex strength

	VortexSearch( vortex, lambda, q_threshold, 1.0 );
	
	PFieldVectorDelete( lambda );

}

/*
 * Searches for the point list on the local domain which contain vortices according to:
 *
 *     \lambda_{ci} >= strength_threshold
 *     | \lambda_{cr} / \lambda_{ci} | <= compactness_threshold
 *
 * Sets the vector vortex_index with the point indices of the global domain
 */
/****************************************************************************************/
void VortexSearch( VortexMap_t &vortex, const PFieldVector_t &lambda, 
		   double strength_threshold, double compactness_threshold )
/****************************************************************************************/
{

	// Clear the in coming vortex
	vortex.clear();

	// Synchronize the needed fields
	lambda[1]->synchronize(); // \lambda_{cr}
	lambda[2]->synchronize(); // \lambda_{ci}

	const size_t N = lambda[0]->getSizeLocal();	

	/// Generate Vortex_t class for the local domain
	Vortex_t *vortex_local = new Vortex_t[N];

#ifdef BOUNDS_CHECK
	// Set the pointer bounds
	vp_lbound = vortex_local;
	vp_ubound = vortex_local + N - 1;
#endif

	//const int *offset_operation = lambda[0]->getOffsetOperation();
	const int *offset_local     = lambda[0]->getOffsetLocal();
	// const int offset_operation_to_global[] = { offset_operation[0] + offset_local[0], 
	// 					   offset_operation[1] + offset_local[1], 
	// 					   offset_operation[2] + offset_local[2] };
#ifdef SANITY_CHECK
	size_t index_check=0;
#endif
	PFIELD_LOOP_LOCAL(lambda[0])

#ifdef SANITY_CHECK
	assert( index_check == _index );
	index_check++;
#endif
	const int ijk_global[] = { _i + offset_local[0], 
				   _j + offset_local[1], 
				   _k + offset_local[2] };

	const size_t index_global = lambda[0]->index( ijk_global[0], ijk_global[1], ijk_global[2] );

#ifdef BOUNDS_CHECK
	assert( 0 <= _index       && _index       < lambda[0]->getSizeLocal() );
	assert( 0 <= index_global && index_global < lambda[0]->getSize() );
#endif

	const double lambda_cr = lambda[1]->getDataLocal()[_index];
	const double lambda_ci = lambda[2]->getDataLocal()[_index];

	// Need to set the volume of the vortex point (cell volume)
	const double cell_volume = lambda[0]->getDXLocal()[ _i ] *
		                   lambda[0]->getDYLocal()[ _j ] *
		                   lambda[0]->getDZLocal()[ _k ];

	//VortexSet( vortex_local[index], index, &(ijk[0]), lambda_ci, lambda_cr / lambda_ci );
	VortexSet( vortex_local[_index], index_global, lambda_ci, lambda_cr / lambda_ci, cell_volume );

	PFIELD_LOOP_END

	// Create the vortex structs for the thresholds
	Vortex_t vortex_strength_threshold; 
        Vortex_t vortex_compactness_threshold_first; 
        Vortex_t vortex_compactness_threshold_last; 

	VortexSet(vortex_strength_threshold,          0, strength_threshold, 0.0,                   0.0);
	VortexSet(vortex_compactness_threshold_first, 0, 0.0,               -compactness_threshold, 0.0);
	VortexSet(vortex_compactness_threshold_last,  0, 0.0,                compactness_threshold, 0.0);

	Vortex_t *vortex_first, *vortex_last;
	size_t N_strong;
	size_t N_vortex;

	// Now sort based on vortex strength
	qsort(vortex_local,N,sizeof(Vortex_t),VortexSortCompareStrength) ;

	// Now have to find the first index that is less than or equal to the threshold strength value
	// Returns the first value (vortex point) that at or above the strength threshold
	Vortex_t *vortex_strong;
	vortex_strong = (Vortex_t *)bsearch( &vortex_strength_threshold, vortex_local+1, 
					     N-1, sizeof(Vortex_t), VortexSearchCompareStrength);

	// Sanity check that the returned vortex is legit
	//assert( vortex_strong != NULL );
	if( vortex_strong == NULL ) 
		goto clean_up;

	// Can now sort based on compactness ratio
	// The number of vortex points that are at or above the strength threshold
	N_strong = (vortex_local + N) - vortex_strong;
	assert( vortex_strong + N_strong == vortex_local + N ); 

	qsort(vortex_strong, N_strong, sizeof(Vortex_t), VortexSortCompareCompactness);

	// If all the vortex points are not compact enough or if they are too
	// compact then there are no vortices that meet the criteria; clean up
	// and exit.
	if( vortex_compactness_threshold_first.compactness > (vortex_strong + N_strong - 1)->compactness ||
	    vortex_compactness_threshold_last.compactness  < vortex_strong->compactness ) {
		printf("no vortices\n");
		goto clean_up;
	} else {

		
		if( vortex_compactness_threshold_first.compactness <= vortex_strong->compactness ) {
			// Take the first vortex point if the threshold
			// compactness is below the compactness of the first
			// vortex point
			vortex_first = vortex_strong; 
		} else {
			vortex_first = (Vortex_t*)bsearch( &vortex_compactness_threshold_first, vortex_strong+1, 
							   N_strong - 1, sizeof(Vortex_t), VortexSearchCompareCompactnessFirst);
			assert(vortex_first != NULL);
		}
	
		if( vortex_compactness_threshold_last.compactness >= (vortex_strong+N_strong-1)->compactness ) {
			// Take the last vortex point if the threshold
			// compactness is above the compactness of the last
			// vortex point
			vortex_last =  vortex_strong + N_strong - 1; 
		} else  {
			vortex_last = (Vortex_t*)bsearch( &vortex_compactness_threshold_last, vortex_strong, 
							  N_strong-1, sizeof(Vortex_t), VortexSearchCompareCompactnessLast);
			assert(vortex_last != NULL);
		}
	}
	
	N_vortex = vortex_last - vortex_first + 1;
	// Finally sort into assending index values
	qsort(vortex_first, N_vortex, sizeof(Vortex_t), VortexSortCompareIndex);

 finalize:
	// Assign the resulting vortex map
	for( size_t n=0; n<N_vortex; n++ ) {
		Vortex_t *v = vortex_first + n;

#ifdef BOUNDS_CHECK
		std::vector<int> _ijk = lambda[0]->ijk(v->index);
		size_t _index_local = lambda[0]->indexLocal( _ijk[0] - lambda[0]->getOffsetLocal()[0], 
							     _ijk[1] - lambda[0]->getOffsetLocal()[1], 
							     _ijk[2] - lambda[0]->getOffsetLocal()[2] );

		assert( 0 <= _index_local && _index_local < lambda[0]->getSizeLocal() );
#endif
		//vortex[v->index] = *v;
		if( !vortex.insert( std::pair<size_t,Vortex_t>(v->index,*v) ).second ) {
			printf("VortexSearch: unable to insert new vortex--already exists\n");
			MPI_Abort(lambda[0]->getMPITopology()->comm,ierr);
		}
	}

 clean_up:

	delete [] vortex_local;
	
}

/****************************************************************************************/
int VortexSortCompareIndex(const void *a,const void *b)
/****************************************************************************************/
{
	if ( ((Vortex_t*)a)->index == ((Vortex_t*)b)->index )
		return 0;
      
	if ( ((Vortex_t*)a)->index < ((Vortex_t*)b)->index ) {
		return -1;
	} else {
		return 1;
	}
}
/****************************************************************************************/
int VortexSortCompareStrength(const void *a,const void *b)
/****************************************************************************************/
{

#ifdef BOUNDS_CHECK
	assert( vp_lbound <= a && a <= vp_ubound );
	assert( vp_lbound <= b && b <= vp_ubound );
#endif

	if ( ((Vortex_t*)a)->strength == ((Vortex_t*)b)->strength )
		return 0;
      
	if ( ((Vortex_t*)a)->strength < ((Vortex_t*)b)->strength ) {
		return -1;
	} else {
		return 1;
	}
}

/****************************************************************************************/
int VortexSortCompareCompactness(const void *a,const void *b)
/****************************************************************************************/
{
	if ( ((Vortex_t*)a)->compactness == ((Vortex_t*)b)->compactness )
		return 0;

	if ( ((Vortex_t*)a)->compactness < ((Vortex_t*)b)->compactness ) {
		return -1;
	} else {
		return 1;
	}
}

/****************************************************************************************/
int  VortexSearchCompareStrength(const void *a,const void *b)
/****************************************************************************************/
{
	if ( ((Vortex_t*)b-1)->strength < ((Vortex_t*)a)->strength && 
	     ((Vortex_t*)a)->strength <= ((Vortex_t*)b)->strength ) 
		return 0;

	if ( ((Vortex_t*)a)->strength < ((Vortex_t*)b-1)->strength ) {
		return -1;
	} else {
		return 1;
	}

}

/****************************************************************************************/
int  VortexSearchCompareCompactnessFirst(const void *a,const void *b)
/****************************************************************************************/
{
	if ( ((Vortex_t*)b-1)->compactness < ((Vortex_t*)a)->compactness && 
	     ((Vortex_t*)a)->compactness <= ((Vortex_t*)b)->compactness ) 
		return 0;

	if ( ((Vortex_t*)a)->compactness < ((Vortex_t*)b-1)->compactness ) {
		return -1;
	} else {
		return 1;
	}

}

/****************************************************************************************/
int  VortexSearchCompareCompactnessLast(const void *a,const void *b)
/****************************************************************************************/
{
	if ( ((Vortex_t*)b)->compactness <= ((Vortex_t*)a)->compactness && 
	     ((Vortex_t*)a)->compactness <  ((Vortex_t*)b+1)->compactness ) 
		return 0;

	if ( ((Vortex_t*)a)->compactness < ((Vortex_t*)b)->compactness ) {
		return -1;
	} else {
		return 1;
	}

}

/****************************************************************************************/
void VortexRegionInit( VortexRegion_t &vortex_region ) 
/****************************************************************************************/
{
	vortex_region.tag = -1;
	vortex_region.vortex_list.clear();
}

/****************************************************************************************/
void VortexRegionSet( VortexRegion_t &vortex_region, const size_t &tag, 
		      VortexMap_t &vortex_list )
/****************************************************************************************/
{
	vortex_region.tag = tag;
	for (VortexMap_t::iterator v = vortex_list.begin(); v != vortex_list.end(); v++ ) 
		vortex_region.vortex_list[v->first] = v->second;
}

// void VortexRegionAppend( VortexRegion_t &vortex_region, const Vortex_t &vortex )
// {
// 	vortex_region.vortex_list.push_back(vortex);
// }

// void VortexRegionCopy( VortexRegion_t &a, const VortexRegion_t &b )
// {
// 	VortexRegionSet( a, b.tag, b.vortex_list );
// }

/****************************************************************************************/
double VortexRegionGetVolume( VortexRegion_t &vortex_region ) 
/****************************************************************************************/
{
	double volume=0.0L;
	for (VortexMap_t::iterator v = vortex_region.vortex_list.begin(); v != vortex_region.vortex_list.end(); v++ ) 
		volume += v->second.volume;

	return volume;
}

/****************************************************************************************/
void VortexRegionComputeVolume( VortexRegion_t &vortex_region ) 
/****************************************************************************************/
{
	vortex_region.volume=VortexRegionGetVolume( vortex_region );
	return;
}

/****************************************************************************************/
void VortexRegionComputeInertia( VortexRegion_t &vortex_region, const PField &host ) 
/****************************************************************************************/
{
	const int *offset_local = host.getOffsetLocal();

	const double *x_local = host.getXLocal();
	const double *y_local = host.getYLocal();
	const double *z_local = host.getZLocal();

	static const int identity[3][3] = { {1, 0, 0}, 
					    {0, 1, 0},
					    {0, 0, 1} };
			

	Inertia_t *inertia = &vortex_region.inertia;

	// Set the inertia tensor to zero
	std::fill_n( *inertia->tensor, 9, 0.0);
	//double inertia_tensor_norm=0.0;

	assert( vortex_region.vortex_list.size() != 0 );

	// Loop through each vortex point in the vortex region
	for (VortexMap_t::iterator v = vortex_region.vortex_list.begin(); v != vortex_region.vortex_list.end(); v++ ) {
		std::vector<int> ijk = host.ijk( v->second.index );
		const int ijk_local[] = { ijk[0] - offset_local[0], 
					  ijk[1] - offset_local[1],
					  ijk[2] - offset_local[2] };
		const double xyz_local[] = { x_local[ ijk_local[0] ], 
					     y_local[ ijk_local[1] ],
					     z_local[ ijk_local[2] ] };
		
		const double r[3] = { xyz_local[0] - vortex_region.barycenter[0],
				      xyz_local[1] - vortex_region.barycenter[1],
				      xyz_local[2] - vortex_region.barycenter[2] };
		const double r2 = pow(r[0],2) + pow(r[1],2) + pow(r[2],2);

		for( int i=0; i<3; i++ ) {
			for( int j=0; j<3; j++ ) {
				inertia->tensor[i][j] += ( r2 * identity[i][j] - r[i]*r[j] ) * v->second.volume;
			}
		}
		if( r2 <= ZERO || v->second.volume <= ZERO ) 
			printf("r2, vortex_volume = %15.7e, %15.7e\n", r2, v->second.volume );

		//inertia_tensor_norm += r2 * v->second.volume;
		
	}

	// for( int i=0; i<3; i++ ) {
	// 	for( int j=0; j<3; j++ ) {
	// 		if( inertia_tensor_norm <= ZERO ) {
	// 			inertia->tensor[i][j] = 0.0;
	// 		} else {
	// 			inertia->tensor[i][j] /= inertia_tensor_norm; //vortex_region.volume;
	// 		}
	// 	}
	// }


	// //printf("identity = \n");
	// printf("tensor = \n");
	// for( int i=0; i<3; i++ ) {
	// 	printf("    ");
	// 	for( int j=0; j<3; j++ ) {
	// 		printf("%15.7e ", inertia->tensor[i][j]);
	// 		//printf("%d ", identity[i][j]);
	// 	}
	// 	printf("\n");
	// }


	// Now get the eigen values of the inertia tensor
	double _eigenvalue_imag[3]; // Imaginary eigenvalues; these should be zero
	lapack_int info = LAPACKE_dgeev(LAPACK_ROW_MAJOR,'N','N', 3, &inertia->tensor[0][0], 3, 
					inertia->eigenvalue, _eigenvalue_imag, NULL, 3, NULL, 3);

	assert( info == 0 );

	//printf("inertia->eigenvalue = %15.7e, %15.7e, %15.7e\n", inertia->eigenvalue[0], inertia->eigenvalue[1], inertia->eigenvalue[2]);
	// Compute the volume of the inertia tensor ellipsoid
	//     x^2 / ( phi sqrt( 1 / I1 ) )^2 + y^2 / ( phi sqrt( 1 / I2 ) )^2 + z^2 / ( phi sqrt( 1 / I3 ) )^2 = 1
	// where phi is determined to give an ellipsoid that has the same inertia and I1, I2, I3 are the eigenvalues of the inertia tensor
	if( inertia->eigenvalue[0] <= ZERO ||
	    inertia->eigenvalue[1] <= ZERO ||
	    inertia->eigenvalue[2] <= ZERO ) {
		inertia->ellipsoid_volume = 0.0;
	} else {

		// Have to find phi
		const double _cache1 = ( inertia->eigenvalue[0] + 
					 inertia->eigenvalue[1] + 
					 inertia->eigenvalue[2] ) *
			               sqrt( inertia->eigenvalue[0] ) * 
			               sqrt( inertia->eigenvalue[1] ) * 
			               sqrt( inertia->eigenvalue[2] ); 
		const double _cache2 = 8.0 / 15.0 * PI * ( 1.0 / inertia->eigenvalue[0] + 
							   1.0 / inertia->eigenvalue[1] +
							   1.0 / inertia->eigenvalue[2] );

		const double phi = pow( _cache1 / _cache2, 1.0/5.0 );

		const double a = phi / sqrt( inertia->eigenvalue[0] );
		const double b = phi / sqrt( inertia->eigenvalue[1] );
		const double c = phi / sqrt( inertia->eigenvalue[2] );

		inertia->ellipsoid_volume = 4.0/3.0 * PI * a * b *c;

	}

	return;
}


/****************************************************************************************/
void VortexRegionComputeBarycenter( VortexRegion_t &vortex_region, const PField &host ) 
/****************************************************************************************/
{
	const int *offset_local = host.getOffsetLocal();

	const double *x_local = host.getXLocal();
	const double *y_local = host.getYLocal();
	const double *z_local = host.getZLocal();

	vortex_region.barycenter.clear();
	vortex_region.barycenter.resize(3,0.0);

	for (VortexMap_t::iterator v = vortex_region.vortex_list.begin(); v != vortex_region.vortex_list.end(); v++ ) {
		std::vector<int> ijk = host.ijk( v->second.index );
		const int ijk_local[] = { ijk[0] - offset_local[0], 
					  ijk[1] - offset_local[1],
					  ijk[2] - offset_local[2] };
		const double xyz_local[] = { x_local[ ijk_local[0] ], 
					     y_local[ ijk_local[1] ],
					     z_local[ ijk_local[2] ] };
		for( int i=0; i<3; i++ ) {
			vortex_region.barycenter[i] += v->second.volume * xyz_local[i];
		}
	}

	for( int i=0; i<3; i++ ) {
		vortex_region.barycenter[i] /= vortex_region.volume;
	}

	return;
}


/****************************************************************************************/
void VortexRegionSearch( VortexRegionMap_t &vortex_region, const VortexMap_t &vortex, 
			 const PField &host)
/****************************************************************************************/
{

	// Clear the in coming vortex region
	vortex_region.clear();

	VortexMap_t vortex_worker = vortex; // Make a copy of the vortex map

	//	printf("VortexRegionSearch: 1\n");

	// Now for each vortex in the vector point list 
	while( !vortex_worker.empty() ) {

		//		printf("VortexRegionSearch: 2\n");

		// Create a new vortex_region
		VortexRegion_t vr;

		VortexRegionInit(vr);
		vr.tag = ++vortex_region_tag; // Starts with 1

		//		printf("VortexRegionSearch: 3\n");

		// Get the vortex iterator for the map
		VortexMap_t::iterator vm = vortex_worker.begin();
	
		// Add the vortex
		size_t index = vm->first;           // global index

		if( !vr.vortex_list.insert( std::pair<size_t,Vortex_t>(vm->first,vm->second)).second ) {
			printf("VortexRegionSearch: unable to insert new vortex point--already exists\n");
			MPI_Abort(host.getMPITopology()->comm,ierr);
		}
		//vr.vortex_list.at(index) = vm->second; // vortex point

		// Erase the vortex that has been added to the region
		vortex_worker.erase(vm);

		//		printf("VortexRegionSearch: 4\n");

		// Populates the vortex list in vr for contigous vortex points
		// starting with the first vortex point added above
		VortexRegionSearchNeighbors( vr, index, vortex_worker, host );
		
		//		printf("VortexRegionSearch: 4.1\n");

		// Insert vortex region into the map
		if( !vortex_region.insert( std::pair<size_t,VortexRegion_t>(vr.tag,vr) ).second ) {
			printf("VortexRegionSearch: unable to insert new vortex region--already exists\n");
			MPI_Abort(host.getMPITopology()->comm,ierr);
		}
		//		printf("VortexRegionSearch: 5\n");
		
	}

	
	//printf("%d: calling VortexRegionSynchronize ...\n", host.getMPITopology()->rank);
	// Now have to synchronize the tags across the processes and renumber
	// them using a global convention
	VortexRegionSynchronize( vortex_region, host );
	
}

/****************************************************************************************/
void VortexRegionSearchNeighbors( VortexRegion_t &vortex_region, 
				  const size_t &index, VortexMap_t &vortex_map, 
				  const PField &host ) 
/****************************************************************************************/
{

	// First check that the vortex_map is not empty; if so simply return
	if( vortex_map.empty() ) return;

	// Get the operation domain ijk values
	std::vector<int> ijk = host.ijk( index ); // global index and ijk values

	const int *offset_local = host.getOffsetLocal();

	// Obtain the local ijk values
	// const int ijk_local[] = { ijk[0] - offset_local[0], 
	//                           ijk[1] - offset_local[1],
	// 			  ijk[2] - offset_local[2] };
	const int *dims_local = host.getDimsLocal(); 

	std::deque<point_t> search_points(1);

	// Current point
	search_points[0].index = index; // global index
	memcpy( search_points[0].ijk, &ijk[0], sizeof(int)*3 ); // global ijk values

	while( search_points.size() > 0 ) {

		// if( host.getMPITopology()->rank == 0 )
		// 	printf("%d: number of search points = %zd for tag %zd\n", 
		// 	       host.getMPITopology()->rank, search_points.size(), vortex_region.tag );

		// new search points
		std::deque<point_t> new_search_points;

		for( std::deque<point_t>::iterator _sp = search_points.begin(); _sp != search_points.end(); _sp++ ) {

			std::deque<point_t> neighbor_points;

			// Need local ijk values to check location in
			// the local domain
			int _ijk_local[] = { _sp->ijk[0] - offset_local[0],
					     _sp->ijk[1] - offset_local[1],
					     _sp->ijk[2] - offset_local[2] };


			point_t _p; // neighbor point to check

			for( int i = _ijk_local[0]-1; i <= _ijk_local[0]+1; i++ ) {
				for( int j = _ijk_local[1]-1; j <= _ijk_local[1]+1; j++ ) {
					for( int k = _ijk_local[2]-1; k <= _ijk_local[2]+1; k++ ) {
						if( host.inDomainLocal( i, j, k ) ) {
							_p.ijk[0] = i + offset_local[0];
							_p.ijk[1] = j + offset_local[1];
							_p.ijk[2] = k + offset_local[2];
							_p.index = host.index( _p.ijk[0], _p.ijk[1], _p.ijk[2] );
							neighbor_points.push_back( _p );
						}
					}
				}
			}

			// if( host.getMPITopology()->rank == 0 )
			// 	printf("%d: number of neighboring points = %zd\n", host.getMPITopology()->rank, neighbor_points.size() );
			// Now check the new neighbor points
			VortexRegionCheckPoints( vortex_region, vortex_map, neighbor_points, host );
			// if( host.getMPITopology()->rank == 0 )
			// 	printf("%d: number of neighboring vortex points = %zd\n", host.getMPITopology()->rank, neighbor_points.size() );

			// Append the new neighbor points from the current sample point
			for( std::deque<point_t>::iterator _np = neighbor_points.begin(); _np != neighbor_points.end(); _np++ )
				new_search_points.push_back( *_np );
			//new_search_points.insert( new_search_points.end(), neighbor_points.begin(), neighbor_points.end() );
		}

		// Now the set the new search points
		search_points.clear();
		search_points = new_search_points;
		new_search_points.clear();

	}

}

void VortexRegionCheckPoints( VortexRegion_t &vortex_region, VortexMap_t &vortex_map, 
				 std::deque<point_t> &points, const PField &host ) {

	std::deque<point_t> vortex_points;
	for( std::deque<point_t>::iterator _p = points.begin(); _p != points.end(); _p++ ) {
				
		VortexMap_t::iterator _vm = vortex_map.find( _p->index ); // Check if the point is a vortex point (global index)
				
		if( _vm != vortex_map.end() ) {
			// Found one; perform the insertion and make sure that it does not already exist in the RB tree
			// printf("%d: VortexRegionSearchPoints: inserting vortex at %zd into region %zd\n",
			//        host.getMPITopology()->rank,_p->index, vortex_region.tag );
			if( !vortex_region.vortex_list.insert( std::pair<size_t,Vortex_t>(_p->index,_vm->second) ).second ) {
				printf("VortexRegionSearchPoints: unable to insert new vortex--already exists\n");
				MPI_Abort( host.getMPITopology()->comm, ierr );
			}
			vortex_points.push_back( *_p ); // Add the point 
			vortex_map.erase(_vm);
			
			if( vortex_map.empty() ) // no more items in the vortex map
				goto finalize;
		}
				
	}

 finalize:
	points.clear();
	points = vortex_points;
			
}


/****************************************************************************************/
void VortexRegionSynchronize( VortexRegionMap_t &vortex_region, const PField &host )
/****************************************************************************************/
{

	const MPITopology_t *mpi_topology = host.getMPITopology();

	MPI_Status status;
	const MPI_Comm comm = mpi_topology->comm;
	const int rank = mpi_topology->rank;

	////////////////////////////////////////////////////////////////////////////////
	//
	// THIS SYNCHRONIZATION PROCEDURE IS BROKEN AND NEEDS
	// FIXING. SET TO ONLY WORK WITH A SINGLE PROCESS FOR NOW
	// 
	// Set the global volume to the local one
	if( mpi_topology->nproc != 1 &&
	    mpi_topology->rank == 0 ) {
		printf("%d: VortexRegionSynchronize: currently only works with single process--setting global volume to local and not synching\n", 
		       mpi_topology->rank);
	}
	for( VortexRegionMap_t::iterator vr=vortex_region.begin(); vr != vortex_region.end(); vr++ ) {
		VortexRegionComputeVolume( vr->second );
		VortexRegionComputeBarycenter( vr->second, host );
		VortexRegionComputeInertia( vr->second, host );
	}

	return;
	////////////////////////////////////////////////////////////////////////////////
	

	//	printf("%d, -1...\n", rank);
	// First create the index to tag map for the vortex region
	std::map<size_t,size_t> vr_index_to_tag; // for a given vortex point with global index, the tag is given
	VortexRegionIndexToTag( vortex_region, vr_index_to_tag ); // Generates mapping with unique tag values

	//printf("%d: 0...\n", rank);

	if( host.getFieldDecomp() != FIELD_DECOMP_SLAB ) {
		std::cout << "VortexRegionSynchronize: currently only FIELD_DECOMP_SLAB decomposition supported\n";
		MPI_Abort(comm,ierr);
	}

	int imin, imax;
	int request=0;
	int tag_global=0;

	const int *dims_operation   = host.getDimsOperation();
	const int *offset_operation = host.getOffsetOperation();
	const int *offset_local     = host.getOffsetLocal();
	const int rind_size         = host.getRindSize();

	size_t N;
	std::map<size_t,size_t> rind_vr_tag_to_index;
	std::vector<size_t>     rind_vr_index_and_tag(1);

	VortexRegionMap_t vortex_region_retag; 

	MPI_Barrier( comm );

	if( mpi_topology->coords[0] != 0 ) {

		////////////////////////////////////////////////////////////////////////////////
		// 1. Receive updated rind region from "prev" process
		//////////////////////////////////////////////////////////////////////////////// 

		//printf("%d: 1...\n", rank);
		// First receive the size of the rind vortex region tag to index vector
		// Wait for the global tag to return
		//printf("%d: receiving N from %d\n", rank, mpi_topology->neighbor_prev[0]);
		MPI_Recv(&N, 1, MPI_UNSIGNED_LONG, mpi_topology->neighbor_prev[0], 1, mpi_topology->comm, &status);
		//printf("%d: received N=%zd from %d\n", rank, N, mpi_topology->neighbor_prev[0]);
		// Resize the rind index and tag vector
		rind_vr_index_and_tag.resize(N);
		// Wait for the vector
 		//printf("%d: receiving rind_vr_index_and_tag from %d\n", rank, mpi_topology->neighbor_prev[0]);
		MPI_Recv(&rind_vr_index_and_tag[0], N, MPI_UNSIGNED_LONG, mpi_topology->neighbor_prev[0], 2, mpi_topology->comm, &status);
 		//printf("%d: received rind_vr_index_and_tag from %d\n", rank, mpi_topology->neighbor_prev[0]);
		////////////////////////////////////////////////////////////////////////////////
		// 2. Update rind region received from "prev"
		//////////////////////////////////////////////////////////////////////////////// 

		//printf("%d: 2...\n", rank);		
		// Now loop through the rind vortex region index and tag vector
		for( size_t n=0; n < N; n+=2 ) {

			size_t index = rind_vr_index_and_tag.at(n);
			size_t tag   = rind_vr_index_and_tag.at(n+1); // New tag

			size_t tag_orig = vr_index_to_tag.at(index);

			// Now modify the original vortex region tag
			VortexRegionMap_t::iterator vr = vortex_region.find(tag_orig); 
			if( vr == vortex_region.end() ) {
				printf("%d: VortexRegionSynchronize: warning vortex region %zd not found in local domain\n",mpi_topology->rank,tag_orig);
				//MPI_Abort(mpi_topology->comm, ierr);
			} else {
				vr->second.tag = tag; // modifying the origin tag value of the
				// vortex region in the origin map with
				// the new tag

				// Update the new map with the updated vortex region
				vortex_region_retag[tag] = vr->second;

				// Now delete the vortex region from the origin vortex region map
				vortex_region.erase(vr);
			}

		}

		////////////////////////////////////////////////////////////////////////////////
		// 3. Retrieve and update global tag index
		//////////////////////////////////////////////////////////////////////////////// 

		//printf("%d: 3...\n", rank);

		// Send request for global tag
		request = 1;
		MPI_Send(&request, 1, MPI_INT, 0, 11, mpi_topology->comm);

		// Wait for the global tag to return
		MPI_Recv(&tag_global, 1, MPI_UNSIGNED_LONG, 0, 12, mpi_topology->comm, &status);
	  
		// size_t tag_global_start = tag_global + 1;

		// // Reserve next set of tags from the remaining number of vortex
		// // regions in the original list;
		// tag_global += vortex_region.size();


		////////////////////////////////////////////////////////////////////////////////
		// 5 Update all other vortex region tags
		//////////////////////////////////////////////////////////////////////////////// 

		//		printf("%d: 5...\n", rank);

		while( !vortex_region.empty() ) {
			
			VortexRegionMap_t::iterator vr = vortex_region.begin(); 
			vr->second.tag = ++tag_global; // modifying the origin tag value of the
			// vortex region in the origin map with
			// the new tag

			// Update the new map with the updated vortex region
			//vortex_region_retag[tag_global_start++] = vr->second;			
			if( !vortex_region_retag.insert( std::pair<size_t,VortexRegion_t>(vr->second.tag,vr->second ) ).second ) {
				printf("%d: VortexRegionSynchronize: unable to insert new vortex region--already exists: dropping vortex\n", 
				       host.getMPITopology()->rank);
				//MPI_Abort(host.getMPITopology()->comm,ierr);
			}
			// Delete the vortex region from the map
			vortex_region.erase(vr);
		}

		// Reassign the vortex region map
		vortex_region = vortex_region_retag;

		vortex_region_retag.clear();

		// Regenerate the vortex region index to tag mapping with the new tags; used for packing 
		vr_index_to_tag.clear(); // clear the old map
		VortexRegionIndexToTag( vortex_region, vr_index_to_tag ); // create the new map for generating the rind_vr_index_to_tag map

		////////////////////////////////////////////////////////////////////////////////
		// 4 Send the updated global tag index back
		//////////////////////////////////////////////////////////////////////////////// 
	 
 		//printf("%d: 4...\n", rank);

		// Send back the global tag
		MPI_Send(&tag_global, 1, MPI_UNSIGNED_LONG, 0, 13, mpi_topology->comm);	  

		// Send request that we are done
		request = 0;
		MPI_Send(&request, 1, MPI_INT, 0, 11, mpi_topology->comm);


	}


	////////////////////////////////////////////////////////////////////////////////
	// 6. Send updated rind region to "next" process
        //////////////////////////////////////////////////////////////////////////////// 

	//	printf("%d: 6...\n", rank);

	// Now loop over the rind region if sending to "next" process
	if( mpi_topology->coords[0] < mpi_topology->dims[0] - 1 ) {

		imin = dims_operation[0] - rind_size;
		imax = imin + rind_size - 1;

		for( int i=imin; i<=imax; i++ ) {
			for (int j=0; j<dims_operation[1]; j++ ) {
				for (int k=0; k<dims_operation[2]; k++ ) {

					int ijk[] = { i + offset_operation[0] + offset_local[0],
						      j + offset_operation[1] + offset_local[1],
						      k + offset_operation[2] + offset_local[2] };

					// global index
					size_t index = host.index( ijk[0], ijk[1], ijk[2] );

					// Now search for the tag in the rind vortex region
					std::map<size_t,size_t>::iterator vr = vr_index_to_tag.find(index);			  
					if( vr != vr_index_to_tag.end() ) {
						// It is a vortex point
						size_t tag = vr->second;
						// Now check if the updated tag has been inserted in to the map
						std::map<size_t,size_t>::iterator rvr = rind_vr_tag_to_index.find(tag);
						if( rvr == rind_vr_tag_to_index.end() ) {
							// The tag,index pair has not been inserted yet
							//rind_vr_tag_to_index[tag]=index;
							if( !rind_vr_tag_to_index.insert( std::pair<size_t,size_t>(tag,index) ).second ) {
								printf("VortexRegionSynchronize: unable to insert new vortex region tag to index map--already exists\n");
								MPI_Abort(host.getMPITopology()->comm,ierr);
							}
						}
					}

				}
			}
		}

		// Create the vector to send to "next" process
		// Now pack the index + tag vector; interleaved; pre-allocated to maximum size
		N = 2*rind_vr_tag_to_index.size();
		rind_vr_index_and_tag.resize(N); // vector

		// Pack the vector to send to the adjacent process
		size_t n=0;
		for( std::map<size_t,size_t>::iterator rvr = rind_vr_tag_to_index.begin(); rvr != rind_vr_tag_to_index.end(); rvr++ ) {
			rind_vr_index_and_tag[n]   = rvr->second; // index
			rind_vr_index_and_tag[n+1] = rvr->first;  // tag
			n+=2;
		};

		// Wait for the global tag to return
		//printf("%d: sending N=%zd to %d\n", rank, N, mpi_topology->neighbor_next[0]);

		// Send the rind index and tag vector
		MPI_Send(&N, 1, MPI_UNSIGNED_LONG, mpi_topology->neighbor_next[0], 1, mpi_topology->comm);
		// Send the rind index and tag vector
		MPI_Send(&rind_vr_index_and_tag[0], N, MPI_UNSIGNED_LONG, mpi_topology->neighbor_next[0], 2, mpi_topology->comm);

	}

	////////////////////////////////////////////////////////////////////////////////
	// Master process serve global tag value to workers
        ////////////////////////////////////////////////////////////////////////////////
	if( mpi_topology->coords[0] == 0 ) {
		// Set the global tag to the number of vortex regions in the
		// master process. This assumes that the votex regions have tags
		// that start at and are incremented by 1
		tag_global = vortex_region.size(); 

		int nworkers = mpi_topology->dims[0] - 1;

		//printf("%d: VortexRegionSynchronize: tag_global = %d, nworkers = %d\n", mpi_topology->rank, tag_global, nworkers);
		
		// Master process
		while( nworkers != 0 ) {
			//printf("%d: VortexRegionSynchronize: nworkers = %d\n", mpi_topology->rank, nworkers);
			// Wait for message
			MPI_Recv(&request, 1, MPI_INT, MPI_ANY_SOURCE, 11, mpi_topology->comm, &status);
			int rank_worker = status.MPI_SOURCE;
			if( request == 1 ) {
				// Send the global tag
				MPI_Send(&tag_global, 1, MPI_UNSIGNED_LONG, rank_worker, 12, mpi_topology->comm);
				// Wait for the global tag to return
				MPI_Recv(&tag_global, 1, MPI_UNSIGNED_LONG, rank_worker, 13, mpi_topology->comm, &status);
			} else {
				// The worker is finished; take it from the pool
				nworkers--;
			}	  
		}

	}

	printf("%d: VortexRegionSynchronize: computing global volume ...\n", mpi_topology->rank);
	// Finally synchronize the vortex volumes
	// Get the total number of vortex regions
	//MPI_Bcast( &tag_global, 1, MPI_UNSIGNED_LONG, 0, mpi_topology->comm );
	size_t N_vr = vortex_region.size();
	size_t N_vr_global;
	MPI_Allreduce( &N_vr, &N_vr_global, 1, MPI_UNSIGNED_LONG, MPI_SUM, mpi_topology->comm);

	// Now everybody creates a vector for each vortex region
	std::vector<double> vr_volume(N_vr_global);
	// Fill the vector with 0.0
	std::fill_n(&vr_volume[0], vr_volume.size(), 0.0L);

	for( VortexRegionMap_t::iterator vr=vortex_region.begin(); vr != vortex_region.end(); vr++ ) {
		const size_t tag = vr->second.tag;
		vr_volume.at(tag-1) = VortexRegionGetVolume( vr->second );
	}

	// Perform the summation of the volumes
	std::vector<double> vr_volume_global( vr_volume.size() );
	MPI_Allreduce ( &vr_volume[0], &vr_volume_global[0], vr_volume.size(), MPI_DOUBLE, MPI_SUM, mpi_topology->comm );

	// Finally set the volume to the global value
	for( VortexRegionMap_t::iterator vr=vortex_region.begin(); vr != vortex_region.end(); vr++ ) {
		const size_t tag = vr->second.tag;
		vr->second.volume = vr_volume_global.at( tag-1 );
	}


	//printf("%d: VortexRegionSynchronize: done\n", mpi_topology->rank);
	
}

/*******************************************************************************/
//
// Generates an (global) index to tag map of the vortex regions.
//
/*******************************************************************************/
/****************************************************************************************/
void VortexRegionIndexToTag( VortexRegionMap_t &vortex_region, 
			     std::map<size_t,size_t> &vr_index_to_tag ) 
/****************************************************************************************/
{

	//printf("WHAT UP\n");
	// Clear the incoming map
	vr_index_to_tag.clear();

	for (VortexRegionMap_t::iterator vr = vortex_region.begin(); vr != vortex_region.end(); vr++ ) {
		for (VortexMap_t::iterator v = vr->second.vortex_list.begin(); v != vr->second.vortex_list.end(); v++ ) {

			size_t index = v->second.index;
			size_t tag   = vr->second.tag;

			if( !vr_index_to_tag.insert( std::pair<size_t,size_t>(index,tag) ).second ) {
				printf("VortexRegionIndexToTag: unable to insert new vortex region index to tag map--already exists\n");
				exit(EXIT_FAILURE);
			}

		}
	}

	//printf("WHAT DOWN\n");

}

/****************************************************************************************/
void VortexRegionTrimBBox( VortexRegionMap_t &_vortex_region, const PField &_host, 
			   const int _irange[2], const int _jrange[2], const int _krange[2] )
/****************************************************************************************/
{

	VortexRegionMap_t vortex_region_trim;
	
	for( VortexRegionMap_t::iterator vr=_vortex_region.begin(); vr != _vortex_region.end(); vr++ ) {
		bool trim_point=false;
		// Loop over the vortex points in the vortex region
		for( VortexMap_t::iterator v=vr->second.vortex_list.begin(); v != vr->second.vortex_list.end(); v++ ) {
			const std::vector<int> ijk = _host.ijk( v->second.index );
			if( ijk[0] < _irange[0] || ijk[0] > _irange[1] || 
			    ijk[1] < _jrange[0] || ijk[1] > _jrange[1] ||
			    ijk[2] < _krange[0] || ijk[2] > _krange[1] ) {
				trim_point = true;
				break;
			}
		}
		if( ! trim_point )
			vortex_region_trim[vr->second.tag] = vr->second;
	}
	_vortex_region.clear();
	_vortex_region = vortex_region_trim;
}

}
