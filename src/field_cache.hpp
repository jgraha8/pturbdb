#ifndef FIELD_CACHE_HPP_
#define FIELD_CACHE_HPP_

#include <vector>
#include <map>
#include <string>

typedef struct {
	std::string name;
	int nfields;
	size_t size;
	std::vector<bool> is_cached;
	std::map<int, double *> data; // Maps the file index to the field
} FieldCache_t;

FieldCache_t *FieldCacheNew( std::string name, int nfields, size_t N );
void FieldCacheDelete(FieldCache_t **field_cache); 
void FieldCacheSet( FieldCache_t *field_cache, std::vector<int> &file_index );

#endif

