#ifndef FIELD_CACHE_HPP_
#define FIELD_CACHE_HPP_

#include <vector>
#include <map>
#include <string>

namespace pturbdb {

template<typename Type> 
class FieldCache
{

private:
	std::string name_;
	int nfields_;
	size_t size_;
	std::vector<bool> is_cached_;
	typename std::map<int, Type *> data_; // Maps the file index to the field

public:	
	FieldCache( std::string name, int nfields, size_t N );
	~FieldCache(); 

	// Getters
	std::string                    &getName()     {return this->name_;}
	int                            &getNFields()  {return this->nfields_;}
	size_t                         &getFieldSize(){return this->field_size_;}
	std::vector<bool>              &getIsCached() {return this->is_cached_;}
	typename std::map<int, Type *> &getData()     {return this->data_;}

	// Setters
	void setCache( std::vector<int> &file_index );

};

template<typename Type> 
FieldCache<Type>::FieldCache( std::string name, int nfields, size_t N )
{
	
	this->name_    = name;
	this->nfields_ = nfields;
	this->size_    = N;
	this->is_cached_.assign(nfields,false);

	for( int i=0; i<nfields; i++ ) {
		// Initialize the data map with bogus keys
		this->data_[-(i+1)] = new Type[N];
	}
}

template<typename Type>
FieldCache<Type>::~FieldCache()
{
	this->is_cached_.clear();
	
	// Clear the map
	typename std::map<int, Type *>::iterator d;
	for( d=this->data_.begin(); d != this->data_.end(); d++ ) {
		delete [] d->second;
	}
	this->data_.clear();
}

/*
 * Sets the data map for the given file index vector. 
 */
template<typename Type>
void FieldCache<Type>::setCache( std::vector<int> &file_index )
{

	if( file_index.size() != (size_t)this->nfields_ ) {
		std::cout << "FieldCacheSet: file index vector not size 4\n";
		exit(EXIT_FAILURE);
	}
	
	// Buffer data vector
	std::vector<Type *> data_buffer(this->nfields_,NULL);
	// Reset the is_cached vector
	this->is_cached_.assign(this->nfields_,false);

	for( int i=0; i<this->nfields_; i++ ) {
		int j=file_index[i];
		typename std::map<int, Type *>::iterator f = this->data_.find(j);
		if( f != this->data_.end() ) {
			// Found the field in the map
			data_buffer[i] = f->second; // Read the value;
			// Drop the entry from the map
			this->data_.erase(f);
			// Set the current field to be cached;
			this->is_cached_[i] = true;
		}
	}
	
	// Make a second sweep updating all NULL data fields;
	for( int i=0; i<this->nfields_; i++ ) {
		if( data_buffer[i] == NULL ) {
			// Assign it to the first entry in the data map
			data_buffer[i] = this->data_.begin()->second;
			// Erase the entry from the map
			this->data_.erase( this->data_.begin() );	
		}
	}

	if( this->data_.size() != 0 ) {
		std::cout << "FieldCacheSet: field extraction error\n";
		exit(EXIT_FAILURE);
	}

	// Finally rebuild the data map
	for( int i=0; i<this->nfields_; i++ ) {
		this->data_[file_index[i]] = data_buffer[i];
	}
	data_buffer.clear();

}

}
#endif

