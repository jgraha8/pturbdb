#ifndef VORTEX_TRACKING_HPP_
#define VORTEX_TRACKING_HPP_

#include <vector>
#include <deque>
#include "core.hpp"

namespace pturbdb {

// Vortex (point) structure
typedef struct {
	size_t index;
	double strength;
	double compactness;
} Vortex_t;

typedef std::map<size_t,Vortex_t> VortexMap_t;

// Vortex region structure
typedef struct {
	size_t tag;
	VortexMap_t vortex_list;
} VortexRegion_t;

typedef std::map<size_t,VortexRegion_t> VortexRegionMap_t;

void VortexInit( Vortex_t &vortex );
void VortexSet( Vortex_t &vortex, size_t index, double strength, double compactness );
void VortexSearch( VortexMap_t &vortex, const PFieldVector_t &lambda, double strength_threshold, double compactness_threshold);
void VortexSearchQ( VortexMap_t &vortex, const PField &Q, double q_threshold );

void VortexRegionInit( VortexRegion_t &vortex_region );
void VortexRegionSet( VortexRegion_t &vortex_region, const size_t &tag, VortexMap_t &vortex_list );
void VortexRegionAppend( VortexRegion_t &vortex_region, const Vortex_t &vortex );
void VortexRegionSearch( VortexRegionMap_t &vortex_region, const VortexMap_t &vortex_list, const PField &host );
void VortexRegionSearchNeighbors( VortexRegion_t &vortex_region, const size_t &vortex_index, 
				  VortexMap_t &vortex_map, const PField &host );
void VortexRegionCheckPoints( VortexRegion_t &vortex_region, VortexMap_t &vortex_map, 
			      std::deque<point_t> &points, const PField &host );
void VortexRegionSynchronize( VortexRegionMap_t &vortex_region, const PField &host );
void VortexRegionIndexToTag( std::map<size_t,size_t> &vr_index_to_tag, VortexRegionMap_t &vortex_region );

int VortexSortCompareIndex(const void *a,const void *b);
int VortexSortCompareStrength(const void *a,const void *b);
int VortexSortCompareCompactness(const void *a,const void *b);

int VortexSearchCompareStrength(const void *a,const void *b);
int VortexSearchCompareCompactnessFirst(const void *a, const void *b);
int VortexSearchCompareCompactnessLast(const void *a, const void *b);

}

#endif
