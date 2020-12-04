
#ifndef RPVG_SRC_PATHCLUSTERS_HPP
#define RPVG_SRC_PATHCLUSTERS_HPP

#include <vector>

#include "sparsepp/spp.h"

#include "paths_index.hpp"
#include "alignment_path.hpp"

using namespace std;


typedef spp::sparse_hash_map<uint32_t, spp::sparse_hash_set<uint32_t> > connected_paths_t;

class PathClusters {

    public: 

        PathClusters(const PathsIndex & paths_index_in, const uint32_t num_threads_in);
        PathClusters(const PathsIndex & paths_index_in, const uint32_t num_threads_in, const spp::sparse_hash_map<vector<AlignmentPath>, uint32_t> & align_paths_index);

        void addReadClusters(const spp::sparse_hash_map<vector<AlignmentPath>, uint32_t> & align_paths_index);
    	void addCallTraversalClusters();

        vector<uint32_t> path_to_cluster_index;
        vector<vector<uint32_t> > cluster_to_paths_index;

        spp::sparse_hash_map<gbwt::SearchState, uint32_t> search_to_cluster_index;

    private: 

        const PathsIndex & paths_index;
        const uint32_t num_threads;

    	void createPathClusters(const connected_paths_t & connected_paths);
        void mergeClusters(const connected_paths_t & connected_clusters);
};

namespace std {

    template<> 
    struct hash<gbwt::SearchState>
    {
        size_t operator()(gbwt::SearchState const & search_state) const
        {
            size_t seed = 0;
            spp::hash_combine(seed, search_state.node);
            spp::hash_combine(seed, search_state.range.first);
            spp::hash_combine(seed, search_state.range.second);

            return seed;
        }
    };
}


#endif

