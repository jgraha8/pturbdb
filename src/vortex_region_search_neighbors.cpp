
// Get the global index
_index = host.index(_i,_j,_k);
// Search for vortex in map
_vm = vortex_map.find( _index );
if( _vm != vortex_map.end() ) {
	// Found one
	if( !vortex_region.vortex_list.insert( std::pair<size_t,Vortex_t>(_index,_vm->second) ).second ) {
		printf("VortexRegionSearchNeighbors: unable to insert new vortex--already exists\n");
		MPI_Abort(host.getMPITopology()->comm,ierr);
	}
	vortex_map.erase(_vm);
	if( vortex_map.empty() ) goto clean_up;
	//VortexRegionSearchNeighbors( vortex_region, _index, vortex_map, host );
 }
