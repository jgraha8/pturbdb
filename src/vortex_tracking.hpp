#ifndef VORTEX_TRACKING_HPP_
#define VORTEX_TRACKING_HPP_

namespace pturbdb {

typedef struct {
	size_t index;
	double strength;
	double compactness;
} Vortex_t;

// public:
// 	Vortex(): index_(0), strength_(0.0), compactness_(0.0) {};
// 	Vortex( size_t index, double strength, double compactness ) 
// 	{
// 		this->index_ = index;
// 		this->strength_ = strength;
// 		this->compactness_ = compactness;
// 	}
// 	~Vortex(){};

// 	size_t getIndex()       const {return this->index_; }
// 	double getStrength()    const {return this->strength_; }
// 	double getCompactness() const {return this->compactness_; }

// 	void setIndex( size_t index ) { this->index_ = index; };
// 	void setStrength( double strength ) { this->strength__ = strength; };
// 	void setCompactness( double compactness ) { this->compactness_ = compactness; };

void VortexSearch( std::vector<size_t> &vortex_index, PFieldVector_t &lambda, PFieldTensor_t &eigvec, double strength_threshold, double compactness_threshold)
int VortexSortCompareStrength(const void *a,const void *b);
int VortexSortCompareCompactness(const void *a,const void *b);

int VortexSearchCompareStrength(const void *a,const void *b);
int VortexSearchCompareCompactnessFirst(const void *a, const void *b);
int VortexSearchCompareCompactnessLast(const void *a, const void *b);

}

#endif
