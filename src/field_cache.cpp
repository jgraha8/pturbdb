#include <string>
#include <map>
#include <vector>
#include <iostream>
#include <cstdlib>

#include "field_cache.hpp"

FieldCache_t *FieldCacheNew( std::string name, int nfields, size_t N )
{
	
	FieldCache_t *field_cache = new FieldCache_t;
	field_cache->name = name;
	field_cache->nfields = nfields;
	field_cache->size = N;
	field_cache->is_cached.assign(4,false);

	for( int i=0; i<nfields; i++ ) {
		// Initialize the data map with bogus keys
		field_cache->data[-(i+1)] = new double[N];
	}
	return field_cache;
}

void FieldCacheDelete(FieldCache_t **field_cache) 
{
	(*field_cache)->is_cached.clear();
	
	// Clear the map
	std::map<int, double *>::iterator d;
	for( d=(*field_cache)->data.begin(); d != (*field_cache)->data.end(); d++ ) {
		delete [] d->second;
	}
	(*field_cache)->data.clear();

	delete *field_cache;
	*field_cache=NULL;
}

/*
 * Sets the data map for the given file index vector. 
 */
void FieldCacheSet( FieldCache_t *field_cache, std::vector<int> &file_index )
{

	if( file_index.size() != (size_t)field_cache->nfields ) {
		std::cout << "FieldCacheSet: file index vector not size 4\n";
		exit(EXIT_FAILURE);
	}
	
	// Buffer data vector
	std::vector<double *> data_buffer(field_cache->nfields,NULL);
	// Reset the is_cached vector
	field_cache->is_cached.assign(field_cache->nfields,false);

	for( int i=0; i<field_cache->nfields; i++ ) {
		int j=file_index[i];
		std::map<int, double *>::iterator f = field_cache->data.find(j);
		if( f != field_cache->data.end() ) {
			// Found the field in the map
			data_buffer[i] = f->second; // Read the value;
			// Drop the entry from the map
			field_cache->data.erase(f);
			// Set the current field to be cached;
			field_cache->is_cached[i] = true;
		}
	}
	
	// Make a second sweep updating all NULL data fields;
	for( int i=0; i<field_cache->nfields; i++ ) {
		if( data_buffer[i] == NULL ) {
			// Assign it to the first entry in the data map
			data_buffer[i] = field_cache->data.begin()->second;
			// Erase the entry from the map
			field_cache->data.erase( field_cache->data.begin() );	
		}
	}

	if( field_cache->data.size() != 0 ) {
		std::cout << "FieldCacheSet: field extraction error\n";
		exit(EXIT_FAILURE);
	}

	// Finally rebuild the data map
	for( int i=0; i<field_cache->nfields; i++ ) {
		field_cache->data[file_index[i]] = data_buffer[i];
	}
	data_buffer.clear();

}
