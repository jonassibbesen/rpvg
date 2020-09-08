
#ifndef RPVG_SRC_PATHCLUSTERS_HPP
#define RPVG_SRC_PATHCLUSTERS_HPP

#include <vector>
#include <unordered_map>
#include <unordered_set>

#include "sparsepp/spp.h"

#include "paths_index.hpp"
#include "alignment_path.hpp"

using namespace std;


typedef spp::sparse_hash_map<uint32_t, spp::sparse_hash_set<uint32_t> > connected_paths_t;

class PathClusters {

    public: 

        PathClusters(const PathsIndex & paths_index_in, const uint32_t num_threads_in);

        void addReadClusters(const spp::sparse_hash_map<vector<AlignmentPath>, uint32_t> & align_paths_index);
    	void addCallTraversalClusters();

        spp::sparse_hash_map<uint32_t, uint32_t> node_to_path_index;

        vector<uint32_t> path_to_cluster_index;
        vector<vector<uint32_t> > cluster_to_paths_index;

    private: 

        const PathsIndex & paths_index;
        const uint32_t num_threads;

    	void createPathClusters(const connected_paths_t & connected_paths);
        void mergeClusters(const connected_paths_t & connected_clusters);

};


#endif

