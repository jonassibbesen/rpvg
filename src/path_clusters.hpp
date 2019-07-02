
#ifndef VGPROB_PATHCLUSTERS_HPP
#define VGPROB_PATHCLUSTERS_HPP

#include <vector>
#include <unordered_map>
#include <unordered_set>

using namespace std;


class PathClusters {

    public: 
    	
    	PathClusters(const unordered_map<int32_t, unordered_set<int32_t> > & connected_paths, const int32_t num_paths);
    	
    	vector<uint32_t> path_to_cluster_index;
    	vector<vector<uint32_t> > cluster_to_path_index;
};


#endif

