
#include "path_clusters.hpp"

#include <assert.h>
#include <queue>
#include <algorithm>


PathClusters::PathClusters(const spp::sparse_hash_map<uint32_t, spp::sparse_hash_set<uint32_t> > & connected_paths, const uint32_t num_paths) {

    path_to_cluster_index = vector<uint32_t>(num_paths, -1);

    for (size_t i = 0; i < num_paths; ++i) {

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

