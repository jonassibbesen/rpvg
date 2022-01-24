
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

        void addNodeClustering(const PathsIndex & paths_index, const bool cluster_strands);
        void createNoteClusterIndex(const PathsIndex & paths_index, const bool cluster_strands);

        vector<uint32_t> path_to_cluster_index;        
        vector<vector<uint32_t> > cluster_to_paths_index;

        vector<uint32_t> path_to_note_cluster_index;

    private: 

        const uint32_t num_threads;
        const uint32_t num_paths;

        void addConnection(vector<spp::sparse_hash_set<uint32_t> > * connections, vector<mutex> * connections_mutexes, const uint32_t anchor_id, const uint32_t id) const;

        vector<spp::sparse_hash_set<uint32_t> > findNodeConnectedPaths(const PathsIndex & paths_index, const bool cluster_strands) const;
        vector<spp::sparse_hash_set<uint32_t> > findNodeConnectedClusters(const PathsIndex & paths_index, const bool cluster_strands) const;

        vector<vector<gbwt::size_type> > getNodePathSets(const size_t node_idx, const PathsIndex & paths_index) const;

    	void createPathClusters(vector<uint32_t> * path_to_cluster, vector<vector<uint32_t> > * cluster_to_paths, const vector<spp::sparse_hash_set<uint32_t> > & connected_paths);
        void mergeClusters(vector<uint32_t> * path_to_cluster, vector<vector<uint32_t> > * cluster_to_paths, const vector<spp::sparse_hash_set<uint32_t> > & connected_clusters);
};


#endif

