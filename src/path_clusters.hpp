
#ifndef RPVG_SRC_PATHCLUSTERS_HPP
#define RPVG_SRC_PATHCLUSTERS_HPP

#include <vector>
#include <unordered_map>
#include <unordered_set>

#include "sparsepp/spp.h"

using namespace std;


class PathClusters {

    public: 
    	
    	PathClusters(const spp::sparse_hash_map<int32_t, spp::sparse_hash_set<int32_t> > & connected_paths, const int32_t num_paths);
    	
    	vector<int32_t> path_to_cluster_index;
    	vector<vector<int32_t> > cluster_to_path_index;
};


#endif
