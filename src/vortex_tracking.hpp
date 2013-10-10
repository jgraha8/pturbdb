#ifndef VORTEX_TRACKING_HPP_
#define VORTEX_TRACKING_HPP_

namespace pturbdb {

typedef struct {
	size_t index;
	double strength;
	double compactness;
} Vortex_t;

void VortexSet( Vortex_t &vortex, size_t index, double strength, double compactness );
void VortexCopy( Vortex_t &a, const Vortex_t &b );
void VortexSearch( std::vector<Vortex_t> &vortex, PFieldVector_t &lambda, PFieldTensor_t &eigvec, double strength_threshold, double compactness_threshold);

int VortexSortCompareIndex(const void *a,const void *b);
int VortexSortCompareStrength(const void *a,const void *b);
int VortexSortCompareCompactness(const void *a,const void *b);

int VortexSearchCompareStrength(const void *a,const void *b);
int VortexSearchCompareCompactnessFirst(const void *a, const void *b);
int VortexSearchCompareCompactnessLast(const void *a, const void *b);

}

#endif
