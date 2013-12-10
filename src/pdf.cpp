
#include <cmath>
#include <cassert>
#include <vector>
#include "pdf.hpp"

void PDFInit( PDF_t &_p, size_t _nbins, double _bin_width, double _abscissa_min ) 
{
	_p.nbins = _nbins;
	_p.bin_width = _bin_width;
	_p.abscissa_min = _abscissa_min;
	PDFResize( _p, _nbins );
}

void PDFInit( PDF_2D_t &_p, const std::vector<size_t> &_nbins, const std::vector<double> &_bin_width, const std::vector<double> &_abscissa_min ) 
{
	assert( _nbins.size() == 2 && _nbins.size() == _bin_width.size() && _nbins.size() == _abscissa_min.size() );

	_p.nbins = _nbins;
	_p.bin_width = _bin_width;
	_p.abscissa_min = _abscissa_min;

	// Intialize other 2 vector variables
	_p.nsamples.resize(2);
	_p.abscissa.resize(2);

	PDFResize( _p, _nbins );
}

void PDFResize( PDF_t &_p, size_t _n ) 
{
	_p.bins.resize( _n + 1, 0 );
	_p.hist.resize( _n + 1, 0.0 );
	_p.pdf.resize( _n + 1, 0.0 );
	_p.abscissa.resize( _n + 1, 0.0 );
	
	_p.nbins = _n + 1;
}

void PDFResize( PDF_2D_t &_p, const std::vector<size_t> &_n ) 
{

	assert( _n.size() == 2 );

	if( _n[0] > 0 ) {
		_p.nbins[0] = _n[0] + 1;
		_p.bins.resize( _n[0]+1 );
		_p.hist.resize( _n[0]+1 );
		_p.pdf.resize( _n[0]+1 );
		_p.abscissa[0].resize( _n[0]+1, 0.0 );
	}

	if( _n[1] > 0 ) {
		_p.nbins[1] = _n[1] + 1;
		for( size_t n=0; n< _p.nbins[0]; n++ ) {
			_p.bins[n].resize( _n[1] + 1, 0 );
			_p.hist[n].resize( _n[1] + 1, 0.0 );
			_p.pdf[n].resize( _n[1] + 1, 0.0 );
		}
		_p.abscissa[1].resize( _n[1] + 1, 0.0 ); 
	}
}


size_t PDFGetSampleBin( PDF_t &_p, double _sample )
{
	return (size_t)floor( ( _sample - _p.abscissa_min ) / _p.bin_width );
}

std::vector<size_t> PDFGetSampleBin( PDF_2D_t &_p, const std::vector<double> &_sample )
{

	assert( _sample.size() == 2 && _sample.size() == _p.bin_width.size() );

	std::vector<size_t> n(2);
	for( int i=0; i<2; i++ ) 
		n[i] = floor( ( _sample[i] - _p.abscissa_min[i] ) / _p.bin_width[i] );

	return n;
}

void PDFBinSample( PDF_t &_p, double _sample ) 
{
	size_t n = PDFGetSampleBin( _p, _sample );

	if( n <  _p.bins.size() ) 
		_p.bins.at(n)++;

};

void PDFBinSample( PDF_2D_t &_p, const std::vector<double> &_sample ) 
{
	
	assert( _sample.size() == 2 );

	std::vector<size_t> n = PDFGetSampleBin( _p, _sample );

	if( n[0] < _p.nbins[0] && n[1] < _p.nbins[1] )
		// Increment the bin
		_p.bins.at(n[0]).at(n[1])++;

};

void PDFComputeNSamples( PDF_t &_p )
{
	_p.nsamples = 0;
	for( std::vector<size_t>::iterator b=_p.bins.begin(); b != _p.bins.end(); b++ ) {
		_p.nsamples += *b;
	}
}


void PDFComputeNSamples( PDF_2D_t &_p )
{
	_p.nsamples_total=0;
	_p.nsamples.clear();

	_p.nsamples.resize( 2 );
	_p.nsamples[0].resize( _p.nbins[0], 0 );
	_p.nsamples[1].resize( _p.nbins[1], 0 );

	for( size_t i=0; i<_p.nbins[0]; i++ ) {
		// _p.nsamples[0][i] = sum(_p.bins[i][:]);
		for( size_t j=0; j<_p.nbins[1]; j++ ) {
			const size_t bins_ = _p.bins[i][j];
			_p.nsamples_total += bins_;
			_p.nsamples[0][i] += bins_;
		}
	}

	for( size_t j=0; j<_p.nbins[1]; j++ ) {
		// _p.nsamples[1][j] = sum(_p.bins[:][j]);
		for( size_t i=0; i<_p.nbins[0]; i++ ) {
			const size_t bins_ = _p.bins[i][j];
			_p.nsamples[1][j] += bins_;
		}
	}

}

void PDFComputeHist( PDF_t &_p ) 
{
	for( size_t i=0; i<_p.nbins; i++ ) {
		_p.hist[i] = _p.bins[i] / _p.bin_width;
	}
};

void PDFComputeHist( PDF_2D_t &_p ) 
{
	const double bw = _p.bin_width[0] * _p.bin_width[1];

	for( size_t i=0; i<_p.nbins[0]; i++ ) {
		for( size_t j=0; j<_p.nbins[1]; j++ ) {
			_p.hist[i][j] = _p.bins[i][j] / bw;
		}
	}
};

void PDFComputePDF( PDF_t &_p ) 
{
	for( size_t i=0; i<_p.nbins; i++ ) {
		_p.pdf[i] = _p.hist[i] / _p.nsamples;
	}
};

void PDFComputePDF( PDF_2D_t &_p ) 
{
	for( size_t i=0; i<_p.nbins[0]; i++ ) {
		for( size_t j=0; j<_p.nbins[1]; j++ ) {
			_p.pdf[i][j] = _p.hist[i][j] / _p.nsamples_total;
		}
	}
};

void PDFComputeAbscissa( PDF_t &_p ) 
{
	size_t n=0;
	for( std::vector<double>::iterator a=_p.abscissa.begin(); a != _p.abscissa.end(); a++ ) {
		*a = ((double)n + 0.5L) * _p.bin_width + _p.abscissa_min;
		n++;
	}
}

void PDFComputeAbscissa( PDF_2D_t &_p ) 
{

	for( int n=0; n<2; n++ ) {
		for( size_t i=0; i<_p.nbins[n]; i++ ) {
			_p.abscissa[n][i] = ((double)i + 0.5L) * _p.bin_width[n] + _p.abscissa_min[n];
		}
	}
}


void PDFComputeAll( PDF_t &_p )
{
	PDFComputeNSamples( _p );
	PDFComputeHist( _p );
	PDFComputePDF( _p );
	PDFComputeAbscissa( _p );
}

void PDFComputeAll( PDF_2D_t &_p )
{
	PDFComputeNSamples( _p );
	PDFComputeHist( _p );
	PDFComputePDF( _p );
	PDFComputeAbscissa( _p );
}
