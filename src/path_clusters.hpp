
#ifndef RPVG_SRC_PATHCLUSTERS_HPP
#define RPVG_SRC_PATHCLUSTERS_HPP

#include <vector>
#include <unordered_map>
#include <unordered_set>

#include "sparsepp/spp.h"

#include "paths_index.hpp"

using namespace std;


class PathClusters {

    public: 

        PathClusters(const uint32_t num_threads_in);

    	void findPathClusters(spp::sparse_hash_map<uint32_t, spp::sparse_hash_set<uint32_t> > * connected_paths, const PathsIndex & paths_index);   
    	void findPathClusters(vector<spp::sparse_hash_map<uint32_t, spp::sparse_hash_set<uint32_t> > > * connected_paths, const PathsIndex & paths_index);

        spp::sparse_hash_map<uint32_t, uint32_t> findPathNodeClusters(const PathsIndex & paths_index);
    	void findCallTraversalClusters(const PathsIndex & paths_index);

        vector<uint32_t> path_to_cluster_index;
        vector<vector<uint32_t> > cluster_to_paths_index;

    private: 

        const uint32_t num_threads;

    	void createPathClusters(const spp::sparse_hash_map<uint32_t, spp::sparse_hash_set<uint32_t> > & connected_paths, const uint32_t num_paths);
};


#endif

