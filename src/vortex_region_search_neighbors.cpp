// Search for vortex in map
_vm = vortex_map.find( _index );
if( _vm != vortex_map.end() ) {
	// Found one
	vortex_region.vortex_list[_index] = _vm->second;
	vortex_map.erase(_vm);
	VortexRegionSearchNeighbors( vortex_region, _index, vortex_map, host );
 }
