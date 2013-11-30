#include <cmath>
#include <cassert>
#include <cstdlib>
#include <cstdio>
#include <cstring>
#include <vector>
#include <deque>
#include <map>
#include "core.hpp"
#include "pfield.hpp"
#include "pfield_math.hpp"
#include "vortex_tracking.hpp"

namespace pturbdb {

size_t vortex_region_tag=0;
	
void VortexInit( Vortex_t &vortex ) {
	vortex.index=-1;
	vortex.strength=0.0L;
	vortex.compactness=0.0L;
}

void VortexSet( Vortex_t &vortex, size_t index, double strength, double compactness ) 
{
	vortex.index = index;
	vortex.strength = strength;
	vortex.compactness = compactness;
}
// void VortexCopy( Vortex_t &a, const Vortex_t &b )
// {
// 	VortexSet( a, b.index, b.strength, b.compactness );
// }

/*
 * Searches for the point list on the local domain which contain vortices according to:
 *
 *     \lambda_{ci} >= strength_threshold
 *     | \lambda_{cr} / \lambda_{ci} | <= compactness_threshold
 *
 * Sets the vector vortex_index with the point indices of the global domain
 */
void VortexSearch( VortexMap_t &vortex, const PFieldVector_t &lambda, const PFieldTensor_t &eigvec, 
		   double strength_threshold, double compactness_threshold )
{

	// Clear the in coming vortex
	vortex.clear();

	// Synchronize the needed fields
	lambda[1]->synchronize(); // \lambda_{cr}
	lambda[2]->synchronize(); // \lambda_{ci}

	const size_t N = lambda[0]->getSizeLocal();	

	/// Generate Vortex_t class for the local domain
	Vortex_t *vortex_local = new Vortex_t[N];

	const int *offset_operation = lambda[0]->getOffsetOperation();
	const int *offset_local     = lambda[0]->getOffsetLocal();
	const int offset_operation_to_global[] = { offset_operation[0] + offset_local[0], 
						   offset_operation[1] + offset_local[1], 
						   offset_operation[2] + offset_local[2] };

	PFIELD_LOOP_OPERATION_TO_LOCAL(lambda[0])

	const int ijk_global[] = { _i + offset_operation_to_global[0], 
				   _j + offset_operation_to_global[1], 
				   _k + offset_operation_to_global[2] };

	const size_t index_global = lambda[0]->index( ijk_global[0], ijk_global[1], ijk_global[2] );

	const double lambda_cr = lambda[1]->getDataLocal()[_index];
	const double lambda_ci = lambda[2]->getDataLocal()[_index];

	//VortexSet( vortex_local[index], index, &(ijk[0]), lambda_ci, lambda_cr / lambda_ci );
	VortexSet( vortex_local[_index], index_global, lambda_ci, lambda_cr / lambda_ci );
	PFIELD_LOOP_END

	// Create the vortex classes for the thresholds
	Vortex_t vortex_strength_threshold; 
        Vortex_t vortex_compactness_threshold_first; 
        Vortex_t vortex_compactness_threshold_last; 

	VortexSet(vortex_strength_threshold, 0, strength_threshold, 0.0);
	VortexSet(vortex_compactness_threshold_first, 0, 0.0, -compactness_threshold);
	VortexSet(vortex_compactness_threshold_last, 0, 0.0, compactness_threshold);

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
			vortex_first = (Vortex_t*)bsearch( &vortex_compactness_threshold_first, vortex_strong+1, 
							   N_strong - 1, sizeof(Vortex_t), VortexSearchCompareCompactnessFirst);
			assert(vortex_first != NULL);
		}
	
		if( vortex_compactness_threshold_last.compactness >= (vortex_strong+N_strong-1)->compactness ) {
			vortex_last =  vortex_strong + N_strong - 1;
		} else  {
			vortex_last = (Vortex_t*)bsearch( &vortex_compactness_threshold_last, vortex_strong, 
							  N_strong-1, sizeof(Vortex_t), VortexSearchCompareCompactnessLast);
			assert(vortex_first != NULL);
		}
	}
	
	// Peform linear search for points with abs(c2) <= threshold.c2 
	N_vortex = vortex_last - vortex_first + 1;
	// Finally sort into assending index values
	qsort(vortex_first, N_vortex, sizeof(Vortex_t), VortexSortCompareIndex);

	for( Vortex_t *v = vortex_first; v != vortex_last + 1; v++ )
		vortex[v->index] = *v;

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


void VortexRegionInit( VortexRegion_t &vortex_region ) {
	vortex_region.tag = -1;
	vortex_region.vortex_list.clear();
}

void VortexRegionSet( VortexRegion_t &vortex_region, const size_t &tag, VortexMap_t &vortex_list )
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

void VortexRegionSearch( VortexRegionMap_t &vortex_region, const VortexMap_t &vortex, const PField &host)
{

	// Clear the in coming vortex region
	vortex_region.clear();

	VortexMap_t vortex_worker = vortex; // Make a copy of the vortex map

	// Now for each vortex in the vector point list 
	while( vortex_worker.size() != 0 ) {

		// Create a new vortex_region
		VortexRegion_t vr;
		VortexRegionInit(vr);
		vr.tag = ++vortex_region_tag; // Starts with 1

		// Get the vortex iterator for the map
		VortexMap_t::iterator vm = vortex_worker.begin();
		
		// Add the vortex
		size_t index = vm->first;
		vr.vortex_list[index] = vm->second;

		// Erase the vortex that has been added to the region
		vortex_worker.erase(vm);

		// Populates the vortex list in vr for contigous vortex points
		// starting with the first vortex point added above
		VortexRegionSearchNeighbors( vr, index, vortex_worker, host );
		
		// Insert vortex region into the map
		vortex_region[vr.tag] = vr;
		
	}

	// Now have to synchronize the tags across the processes and renumber
	// them using a global convention
	VortexRegionSynchronize( vortex_region, host );
	
}

void VortexRegionSearchNeighbors( VortexRegion_t &vortex_region, const size_t &vortex_index, 
				  VortexMap_t &vortex_map, const PField &host ) 
{

	// Get the operation domain ijk values
	const std::vector<int> ijk_global = host.ijk( vortex_index );

	const int *offset_local = host.getOffsetLocal();

	const int ijk_local[] = { ijk_global[0] - offset_local[0],
	                          ijk_global[1] - offset_local[1],
				  ijk_global[2] - offset_local[2] };

	const int *dims_local = host.getDimsLocal();

	size_t _index;
	VortexMap_t::iterator vm;

	// Have to check that the host field has a rind
	if( ijk_local[0] + 1 < dims_local[0] ) {
		// Get the global index
		_index = host.index( ijk_local[0] + 1 + offset_local[0], 
				     ijk_local[1]     + offset_local[1],
				     ijk_local[2]     + offset_local[2] );
#include "vortex_region_search_neighbors.cpp"
	}

	if( 0 <= ijk_local[0] - 1 ) {
		// Get the global index
		_index = host.index( ijk_local[0] - 1 + offset_local[0], 
				     ijk_local[1]     + offset_local[1],
				     ijk_local[2]     + offset_local[2] );
#include "vortex_region_search_neighbors.cpp"
	}

	if( ijk_local[1] + 1 < dims_local[1] ) {
		// Get the global index
		_index = host.index( ijk_local[0]     + offset_local[0], 
				     ijk_local[1] + 1 + offset_local[1],
				     ijk_local[2]     + offset_local[2] );
#include "vortex_region_search_neighbors.cpp"
	}

	if( 0 <= ijk_local[1] - 1 ) {
		// Get the global index
		_index = host.index( ijk_local[0]     + offset_local[0], 
				     ijk_local[1] - 1 + offset_local[1],
				     ijk_local[2]     + offset_local[2] );
#include "vortex_region_search_neighbors.cpp"
	}

	if( ijk_local[2] + 1 < dims_local[2] ) {
		// Get the global index
		_index = host.index( ijk_local[0]     + offset_local[0], 
				     ijk_local[1]     + offset_local[1],
				     ijk_local[2] + 1 + offset_local[2] );
#include "vortex_region_search_neighbors.cpp"
	}

	if( 0 <= ijk_local[2] - 1 ) {

		// Get the global index
		_index = host.index( ijk_local[0]     + offset_local[0], 
				     ijk_local[1]     + offset_local[1],
				     ijk_local[2] - 1 + offset_local[2] );
#include "vortex_region_search_neighbors.cpp"
	}

}

void VortexRegionSynchronize( VortexRegionMap_t &vortex_region, const PField &host )
{
	
}


}
