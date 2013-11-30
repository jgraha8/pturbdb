// Search for vortex in map
vm = vortex_map.find( _index );
if( vm != vortex_map.end() ) {
	// Found one
	vortex_region.vortex_list[_index] = vm->second;
	vortex_map.erase(vm);
	VortexRegionSearchNeighbors( vortex_region, _index, vortex_map, host );
 }
