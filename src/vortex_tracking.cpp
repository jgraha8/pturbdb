#include <cmath>
#include <cassert>
#include <cstdlib>
#include <cstdio>
#include "core.hpp"
#include "pfield.hpp"
#include "pfield_math.hpp"
#include "vortex_tracking.hpp"

namespace pturbdb {

/*
 * Searches for the point list on the local domain which contain vortices according to:
 *
 *     \lambda_{ci} >= strength_threshold
 *     | \lambda_{cr} / \lambda_{ci} | <= compactness_threshold
 *
 * Sets the vector vortex_index with the point indices of the local domain
 */
void VortexSearch( std::vector<size_t> &vortex_index, PFieldVector_t &lambda, PFieldTensor_t &eigvec, double strength_threshold, double compactness_threshold )
{

	// Clear the in comming vortex
	vortex_index.clear();

	// Synchronize the needed fields
	lambda[1]->synchronize(); // \lambda_{cr}
	lambda[2]->synchronize(); // \lambda_{ci}

	const size_t N = lambda[0]->getSizeLocal();

	/// Generate Vortex_t class for the local domain
	Vortex_t *vortex = new Vortex_t[N];

	for( size_t n=0; n<N; n++ ) {
		double _cache = lambda[2]->getDataLocal()[n];
		vortex[n] = Vortex_t(n, _cache, lambda[1]->getDataLocal()[n] / _cache );
	}

	// Create the vortex classes for the thresholds
	Vortex_t vortex_strength_threshold(0, strength_threshold, 0.0);
        Vortex_t vortex_compactness_threshold_first(0,0.0,-compactness_threshold);
        Vortex_t vortex_compactness_threshold_last(0,0.0,compactness_threshold);
	Vortex_t *vortex_first, *vortex_last;

	// Now sort based on vortex strength
	qsort(vortex,N,sizeof(Vortex_t),Vortex_tSortCompareStrength) ;

	// Now have to find the first index that is less than or equal to the threshold strength value
	// Returns the first value that at or above the threshold
	Vortex_t *vortex_strong = (Vortex_t *)bsearch( &vortex_strength_threshold, vortex+1, N-1, sizeof(Vortex_t), VortexSearchCompareStrength);

	if( vortex_strong == NULL ) goto clean_up;

	// Can now sort based on compactness ratio
	const size_t N_strong = (vortex + N) - vortex_strong;
	assert( vortex_strong + N_strong - 1 == vortex + N - 1 );
	qsort(vortex_strong, N_strong, sizeof(Vortex_t), VortexSortCompareCompactness);


	if( vortex_compactness_threshold_first.compactness > (vortex_strong + N_strong - 1)->compactness &&
	    vortex_compactness_threshold_last.compactness < vortex_strong->compactness ) {
		printf("no vortices\n");
		goto clean_up;
	} else {
		if( vortex_compactness_threshold_first.compactness <= vortex_strong->compactness ) {
			vortex_first = vortex_strong;
		} else {
			vortex_first = (Vortex_t*)bsearch( &vortex_compactness_threshold_first, vortex_strong+1, N_strong - 1, sizeof(Vortex_t), VortexSearchCompareCompactnessFirst);
		}
	
		if( vortex_compactness_threshold_last.compactness >= (vortex_strong+N_strong-1)->compactness ) {
			vortex_last =  vortex_strong + N_strong - 1;
		} else  {
			vortex_last = (Vortex_t*)bsearch( &vortex_compactness_threshold_last, vortex_strong, N_strong-1, sizeof(Vortex_t), VortexSearchCompareCompactnessLast);
		}
	}
	
	// Peform linear search for points with abs(c2) <= threshold.c2 
	const size_t N_vortex = vortex_last - vortex_first + 1;

	// Finally sort into assending index values
	qsort(vortex_first, N_vortex, sizeof(Vortex_t), VortexSortCompareIndex);

	vortex_index.assign(N_vortex,-1);
	for( size_t n=0; n<N_vortex; n++ ) 
		vortex_index[n] = (vortex_first+n)->index; 
		
 clean_up:

	delete [] vortex;
	
}

int VortexSortCompareIndex(const void *a,const void *b)
{
	if ( ((Vortex_t*)a)->index == ((Vortex_t*)b)->index )
		return 0;
      
	if ( ((Vortex_t*)a)->index < ((Vortex_t*)b)->index ) {
		return -1;
	} else {
		return 1;
	}
}
int VortexSortCompareStrength(const void *a,const void *b)
{
	if ( ((Vortex_t*)a)->strength == ((Vortex_t*)b)->strength )
		return 0;
      
	if ( ((Vortex_t*)a)->strength < ((Vortex_t*)b)->strength ) {
		return -1;
	} else {
		return 1;
	}
}

int VortexSortCompareCompactness(const void *a,const void *b)
{
	if ( ((Vortex_t*)a)->compactness == ((Vortex_t*)b)->compactness )
		return 0;

	if ( ((Vortex_t*)a)->compactness < ((Vortex_t*)b)->compactness ) {
		return -1;
	} else {
		return 1;
	}
}

int  VortexSearchCompareStrength(const void *a,const void *b)
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
int  VortexSearchCompareCompactnessFirst(const void *a,const void *b)
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
int  VortexSearchCompareCompactnessLast(const void *a,const void *b)
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

}
