#ifndef PDF_HPP_
#define PDF_HPP_

#include <vector>

typedef struct {
	size_t nbins;
	double bin_width;
	double abscissa_zero;

	size_t nsamples;

	std::vector<size_t> bins;
	std::vector<double> hist;
	std::vector<double> pdf;
	std::vector<double> abscissa;
} PDF_t;

typedef struct {

	std::vector<size_t> nbins;
	std::vector<double> bin_width;
	std::vector<double> abscissa_zero;

	size_t nsamples_total;
	std::vector<std::vector<size_t> > nsamples; 

	std::vector<std::vector<size_t> > bins;
	std::vector<std::vector<double> > hist;
	std::vector<std::vector<double> > pdf;
	std::vector<std::vector<double> > abscissa;

} PDF_2D_t;

void PDFInit( PDF_t &_p, size_t _nbins, double _bin_width, double _abscissa_zero );
void PDFInit( PDF_2D_t &_p, const std::vector<size_t> &_nbins, const std::vector<double> &_bin_width, 
	      const std::vector<double> &_abscissa_zero );

void PDFResize( PDF_t &_p, size_t _n );
void PDFResize( PDF_2D_t &_p, const std::vector<size_t> &_n );

size_t PDFGetSampleBin( PDF_t &_p, double _sample );
std::vector<size_t> PDFGetSampleBin( PDF_2D_t &_p, const std::vector<double> &_sample );

void PDFBinSample( PDF_t &_p, double _sample );
void PDFBinSample( PDF_2D_t &_p, const std::vector<double> &_sample );

void PDFComputeNSamples( PDF_t &_p );
void PDFComputeNSamples( PDF_2D_t &_p );

void PDFComputeHist( PDF_t &_p ) ;
void PDFComputeHist( PDF_2D_t &_p ) ;

void PDFComputePDF( PDF_t &_p ) ;
void PDFComputePDF( PDF_2D_t &_p ) ;

void PDFComputeAbscissa( PDF_t &_p ) ;
void PDFComputeAbscissa( PDF_2D_t &_p ) ;

void PDFComputeAll( PDF_t &_p );
void PDFComputeAll( PDF_2D_t &_p );

#endif
