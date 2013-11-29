#include <cmath>
#include <cassert>
#include <cstdlib>
#include <cstdio>
#include <vector>
#include <deque>
#include <map>
#include "core.hpp"
#include "pfield.hpp"
#include "pfield_math.hpp"
#include "vortex_tracking.hpp"

namespace pturbdb {


void VortexSet( Vortex_t &vortex, size_t index, double strength, double compactness ) 
{
	vortex.index = index;
	vortex.strength = strength;
	vortex.compactness = compactness;
}
void VortexCopy( Vortex_t &a, const Vortex_t &b )
{
	VortexSet( a, b.index, b.strength, b.compactness );
}

/*
 * Searches for the point list on the local domain which contain vortices according to:
 *
 *     \lambda_{ci} >= strength_threshold
 *     | \lambda_{cr} / \lambda_{ci} | <= compactness_threshold
 *
 * Sets the vector vortex_index with the point indices of the local domain
 */
void VortexSearch( std::vector<Vortex_t> &vortex, const PFieldVector_t &lambda, const PFieldTensor_t &eigvec, 
		   double strength_threshold, double compactness_threshold )
{

	// Clear the in coming vortex
	vortex.clear();

	// Synchronize the needed fields
	lambda[1]->synchronize(); // \lambda_{cr}
	lambda[2]->synchronize(); // \lambda_{ci}

	//	const size_t N = lambda[0]->getSizeLocal();	
	const size_t N = lambda[0]->getSizeOperation();

	/// Generate Vortex_t class for the local domain
	Vortex_t *vortex_local = new Vortex_t[N];

	PFIELD_LOOP_OPERATION_TO_LOCAL(lambda[0])
	//for( size_t n=0; n<N; n++ ) {
		// double _cache = lambda[2]->getDataLocal()[n];
		// VortexSet( vortex_local[n], n, _cache, lambda[1]->getDataLocal()[n] / _cache );
        double _cache = lambda[2]->getDataLocal()[_index];
	size_t _index_operation = lambda[0]->indexOperation(_i,_j,_k);
	VortexSet( vortex_local[_index_operation], _index_operation, _cache, lambda[1]->getDataLocal()[_index] / _cache );
		//	}
        PFIELD_LOOP_END

	// Create the vortex classes for the thresholds
	Vortex_t vortex_strength_threshold; VortexSet(vortex_strength_threshold, 0, strength_threshold, 0.0);
        Vortex_t vortex_compactness_threshold_first; VortexSet(vortex_compactness_threshold_first, 0, 0.0, -compactness_threshold);
        Vortex_t vortex_compactness_threshold_last; VortexSet(vortex_compactness_threshold_last, 0, 0.0, compactness_threshold);
	Vortex_t *vortex_first, *vortex_last;
	size_t N_strong;
	size_t N_vortex;

	// Now sort based on vortex strength
	qsort(vortex_local,N,sizeof(Vortex_t),VortexSortCompareStrength) ;

	// Now have to find the first index that is less than or equal to the threshold strength value
	// Returns the first value that at or above the threshold
	Vortex_t *vortex_strong = (Vortex_t *)bsearch( &vortex_strength_threshold, vortex_local+1, N-1, sizeof(Vortex_t), VortexSearchCompareStrength);

	assert( vortex_strong != NULL );

	// Can now sort based on compactness ratio
	N_strong = (vortex_local + N) - vortex_strong;
	assert( vortex_strong + N_strong - 1 == vortex_local + N - 1 );
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
			assert(vortex_first != NULL);
		}
	
		if( vortex_compactness_threshold_last.compactness >= (vortex_strong+N_strong-1)->compactness ) {
			vortex_last =  vortex_strong + N_strong - 1;
		} else  {
			vortex_last = (Vortex_t*)bsearch( &vortex_compactness_threshold_last, vortex_strong, N_strong-1, sizeof(Vortex_t), VortexSearchCompareCompactnessLast);
			assert(vortex_first != NULL);
		}
	}
	
	// Peform linear search for points with abs(c2) <= threshold.c2 
	N_vortex = vortex_last - vortex_first + 1;
	// Finally sort into assending index values
	qsort(vortex_first, N_vortex, sizeof(Vortex_t), VortexSortCompareIndex);

	vortex.resize(N_vortex);
	for( size_t n=0; n<N_vortex; n++ ) 
		VortexCopy( vortex[n], *(vortex_first+n) );
		
 clean_up:

	delete [] vortex_local;
	
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

void VortexRegionSet( VortexRegion_t &vortex_region, const size_t &tag, const std::deque<Vortex_t> &vortex_list )
{
	vortex_region.tag = tag;
	vortex_region.vortex_list = vortex_list;
}

void VortexRegionAppend( VortexRegion_t &vortex_region, const Vortex_t &vortex )
{
	vortex_region.vortex_list.push_back(vortex);
}

void VortexRegionCopy( VortexRegion_t &a, const VortexRegion_t &b )
{
	a.tag = b.tag;
	a.vortex_list = b.vortex_list;
}

void VortexRegionSearch( std::vector<VortexRegion_t> &vortex_region, const std::vector<Vortex_t> &vortex )
{

	// Clear the in coming vortex region
	vortex_region.clear();

	// First generate a map of the vortex points Maps the
	// operation grid index to the vortex index such that for a
	// given operation grid index, the vortex at the index can be
	// retrieved with: 
	//     vortex[vortex_map[index]]
	std::map<size_t,size_t> vortex_map; 

	for( size_t i=0; i<vortex.size(); i++ ) 
		vortex_map[vortex[i].index]=i;
    


}


}
