
#ifndef RPVG_SRC_PATHCLUSTERS_HPP
#define RPVG_SRC_PATHCLUSTERS_HPP

#include <vector>

#include "sparsepp/spp.h"

#include "paths_index.hpp"
#include "alignment_path.hpp"

using namespace std;


class PathClusters {

    public: 

        PathClusters(const uint32_t num_threads_in, const PathsIndex & paths_index, const spp::sparse_hash_map<vector<AlignmentPath>, uint32_t> & align_paths_index);

        void addNodeClusters(const PathsIndex & paths_index);

        vector<uint32_t> path_to_cluster_index;
        vector<vector<uint32_t> > cluster_to_paths_index;

    private: 

        const uint32_t num_threads;
        const uint32_t num_paths;

    	void createPathClusters(const vector<spp::sparse_hash_set<uint32_t> > & connected_paths);
        void mergeClusters(const vector<spp::sparse_hash_set<uint32_t> > & connected_clusters);
};


#endif

