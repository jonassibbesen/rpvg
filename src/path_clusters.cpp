
#include <assert.h>
#include <queue>
#include <algorithm>

#include "path_clusters.hpp"


PathClusters::PathClusters(const spp::sparse_hash_map<uint32_t, spp::sparse_hash_set<uint32_t> > & connected_paths, const uint32_t num_paths) {

    findPathClusters(connected_paths, num_paths);
}

PathClusters::PathClusters(const vector<spp::sparse_hash_map<uint32_t, spp::sparse_hash_set<uint32_t> > > & connected_paths, const uint32_t num_paths) {

    auto merged_connected_paths = connected_paths.front();

    for (size_t i = 1; i < connected_paths.size(); ++i) {

        for (auto & path_component: connected_paths.at(i)) {

            auto merged_connected_paths_it = merged_connected_paths.emplace(path_component.first, spp::sparse_hash_set<uint32_t>());

            for (auto & path_id: path_component.second) {

                merged_connected_paths_it.first->second.emplace(path_id);
            }
        }
    }

    findPathClusters(merged_connected_paths, num_paths);
}

void PathClusters::findPathClusters(const spp::sparse_hash_map<uint32_t, spp::sparse_hash_set<uint32_t> > & connected_paths, const uint32_t num_paths) {

    path_to_cluster_index = vector<uint32_t>(num_paths, -1);

    for (uint32_t i = 0; i < num_paths; ++i) {

        if (path_to_cluster_index.at(i) == -1) {

            queue<uint32_t> search_queue;
            search_queue.push(i);

            cluster_to_path_index.emplace_back(vector<uint32_t>());

            while (!search_queue.empty()) {

                auto cur_path = search_queue.front();

                bool is_first_visit = (path_to_cluster_index.at(cur_path) == -1);
                assert(is_first_visit || path_to_cluster_index.at(cur_path) == cluster_to_path_index.size() - 1);

                path_to_cluster_index.at(cur_path) = cluster_to_path_index.size() - 1;
                
                if (is_first_visit) {

                    cluster_to_path_index.back().emplace_back(cur_path);
                    auto connected_paths_it = connected_paths.find(cur_path);

                    if (connected_paths_it != connected_paths.end()) {

                        for (auto & path: connected_paths_it->second) {

                            if (path_to_cluster_index.at(path) == -1) {

                                search_queue.push(path);
                            }
                        }
                    }
                }

                search_queue.pop();
            }

            sort(cluster_to_path_index.back().begin(), cluster_to_path_index.back().end());
        }
    }
}

