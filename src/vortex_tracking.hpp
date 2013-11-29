#ifndef VORTEX_TRACKING_HPP_
#define VORTEX_TRACKING_HPP_

#include <vector>
#include <deque>

namespace pturbdb {

// Vortex (point) structure
typedef struct {
	size_t index;
	double strength;
	double compactness;
} Vortex_t;

// Vortex region structure
typedef struct {
	size_t tag;
	std::deque<Vortex_t> vortex_list;
} VortexRegion_t;

void VortexSet( Vortex_t &vortex, size_t index, double strength, double compactness );
void VortexCopy( Vortex_t &a, const Vortex_t &b );
void VortexSearch( std::vector<Vortex_t> &vortex, const PFieldVector_t &lambda, const PFieldTensor_t &eigvec, 
		   double strength_threshold, double compactness_threshold);

void VortexRegionSet( VortexRegion_t &vortex_region, const size_t &tag, const std::deque<Vortex_t> &vortex_list );
void VortexRegionAppend( VortexRegion_t &vortex_region, const Vortex_t &vortex );
void VortexRegionCopy( VortexRegion_t &a, const VortexRegion_t &b );
void VortexRegionSearch( std::vector<VortexRegion_t> &vortex_region, const std::vector<Vortex_t> &vortex );

int VortexSortCompareIndex(const void *a,const void *b);
int VortexSortCompareStrength(const void *a,const void *b);
int VortexSortCompareCompactness(const void *a,const void *b);

int VortexSearchCompareStrength(const void *a,const void *b);
int VortexSearchCompareCompactnessFirst(const void *a, const void *b);
int VortexSearchCompareCompactnessLast(const void *a, const void *b);

}

#endif
