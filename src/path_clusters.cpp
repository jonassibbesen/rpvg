#include <assert.h>
#include <queue>
#include <algorithm>

#include "path_clusters.hpp"
#include "utils.hpp"


static const uint32_t paths_per_mutex = 100;
static const uint32_t clusters_per_mutex = 100;

PathClusters::PathClusters(const uint32_t num_threads_in, const PathsIndex & paths_index, const spp::sparse_hash_map<vector<AlignmentPath>, uint32_t> & align_paths_index) : num_threads(num_threads_in), num_paths(paths_index.numberOfPaths()) {

    vector<spp::sparse_hash_set<uint32_t> > connected_paths(num_paths, spp::sparse_hash_set<uint32_t>());
    vector<mutex> connected_paths_mutexes(ceil(num_paths / static_cast<double>(paths_per_mutex)));

    #pragma omp parallel num_threads(num_threads)
    {
        #pragma omp for schedule(static, 1)
        for (size_t i = 0; i < num_threads; ++i) {

            auto align_paths_index_it = align_paths_index.begin();
            advance(align_paths_index_it, i);

            uint32_t cur_align_paths_idx = i; 

            while (cur_align_paths_idx < align_paths_index.size()) {

                assert(align_paths_index_it->first.size() > 1);
                assert(align_paths_index_it->first.back().gbwt_search.first.empty());

                uint32_t anchor_path_id = 0;

                for (size_t j = 0; j < align_paths_index_it->first.size() - 1; ++j) {

                    auto align_path_ids = paths_index.locatePathIds(align_paths_index_it->first.at(j).gbwt_search);
                    assert(!align_path_ids.empty());

                    if (j == 0) {

                        anchor_path_id = align_path_ids.front();
                    }

                    auto anchor_path_mutex_idx = floor(anchor_path_id / static_cast<double>(paths_per_mutex));

                    for (auto & path_id: align_path_ids) {

                        if (anchor_path_id != path_id) {

                            auto path_mutex_idx = floor(path_id / static_cast<double>(paths_per_mutex));

                            if (path_mutex_idx != anchor_path_mutex_idx) {

                                lock(connected_paths_mutexes.at(anchor_path_mutex_idx), connected_paths_mutexes.at(path_mutex_idx));

                            } else {

                                connected_paths_mutexes.at(anchor_path_mutex_idx).lock();
                            }
                                                    
                            if (connected_paths.at(anchor_path_id).emplace(path_id).second) {

                                connected_paths.at(path_id).emplace(anchor_path_id);
                            }

                            connected_paths_mutexes.at(anchor_path_mutex_idx).unlock();

                            if (path_mutex_idx != anchor_path_mutex_idx) {

                                connected_paths_mutexes.at(path_mutex_idx).unlock();
                            }
                        }
                    }
                }

                advance(align_paths_index_it, num_threads);
                cur_align_paths_idx += num_threads;
            }
        }
    }

    createPathClusters(connected_paths);
}

void PathClusters::addNodeClusters(const PathsIndex & paths_index) {

    vector<spp::sparse_hash_set<uint32_t> > connected_clusters(cluster_to_paths_index.size(), spp::sparse_hash_set<uint32_t>());
    vector<mutex> connected_clusters_mutexes(ceil(cluster_to_paths_index.size() / static_cast<double>(clusters_per_mutex)));

    #pragma omp parallel num_threads(num_threads)
    {
        #pragma omp for schedule(static)
        for (size_t i = 1; i <= paths_index.numberOfNodes(); ++i) {

            vector<vector<gbwt::size_type> > node_path_id_sets;

            pair<gbwt::SearchState, gbwt::size_type> gbwt_search;
            paths_index.find(&gbwt_search, gbwt::Node::encode(i, false));

            if (!gbwt_search.first.empty()) {

                node_path_id_sets.emplace_back(paths_index.locatePathIds(gbwt_search));
            }

            if (!paths_index.bidirectional()) {

                paths_index.find(&gbwt_search, gbwt::Node::encode(i, true));

                if (!gbwt_search.first.empty()) {

                    node_path_id_sets.emplace_back(paths_index.locatePathIds(gbwt_search));
                }
            }

            for (auto & node_path_ids: node_path_id_sets) {

                if (!node_path_ids.empty()) {

                    auto anchor_cluster_id = path_to_cluster_index.at(node_path_ids.front());
                    auto anchor_cluster_mutex_idx = floor(anchor_cluster_id / static_cast<double>(clusters_per_mutex));

                    for (auto & path_id: node_path_ids) {

                        auto cluster_id = path_to_cluster_index.at(path_id);

                        if (anchor_cluster_id != cluster_id) {

                            auto cluster_mutex_idx = floor(cluster_id / static_cast<double>(clusters_per_mutex));

                            if (cluster_mutex_idx != anchor_cluster_mutex_idx) {

                                lock(connected_clusters_mutexes.at(anchor_cluster_mutex_idx), connected_clusters_mutexes.at(cluster_mutex_idx));

                            } else {

                                connected_clusters_mutexes.at(anchor_cluster_mutex_idx).lock();
                            }
                                                    
                            if (connected_clusters.at(anchor_cluster_id).emplace(cluster_id).second) {

                                connected_clusters.at(cluster_id).emplace(anchor_cluster_id);
                            }

                            connected_clusters_mutexes.at(anchor_cluster_mutex_idx).unlock();

                            if (cluster_mutex_idx != anchor_cluster_mutex_idx) {

                                connected_clusters_mutexes.at(cluster_mutex_idx).unlock();
                            }
                        }
                    }
                }
            }
        }
    }

    if (!connected_clusters.empty()) {

        mergeClusters(connected_clusters);
    }
}

void PathClusters::createPathClusters(const vector<spp::sparse_hash_set<uint32_t> > & connected_paths) {

    assert(path_to_cluster_index.empty());
    assert(cluster_to_paths_index.empty());

    path_to_cluster_index = vector<uint32_t>(num_paths, -1);

    for (uint32_t i = 0; i < num_paths; ++i) {

        if (path_to_cluster_index.at(i) == -1) {

            queue<uint32_t> search_queue;
            search_queue.push(i);

            cluster_to_paths_index.emplace_back(vector<uint32_t>());

            while (!search_queue.empty()) {

                auto cur_path = search_queue.front();

                bool is_first_visit = (path_to_cluster_index.at(cur_path) == -1);
                assert(is_first_visit || path_to_cluster_index.at(cur_path) == cluster_to_paths_index.size() - 1);

                path_to_cluster_index.at(cur_path) = cluster_to_paths_index.size() - 1;
                
                if (is_first_visit) {

                    cluster_to_paths_index.back().emplace_back(cur_path);

                    for (auto & next_path: connected_paths.at(cur_path)) {

                        if (path_to_cluster_index.at(next_path) == -1) {

                            search_queue.push(next_path);
                        }
                    }
                }

                search_queue.pop();
            }

            sort(cluster_to_paths_index.back().begin(), cluster_to_paths_index.back().end());
        }
    }
}

void PathClusters::mergeClusters(const vector<spp::sparse_hash_set<uint32_t> > & connected_clusters) {

    auto old_cluster_to_paths_index = cluster_to_paths_index;
    cluster_to_paths_index.clear();

    vector<bool> visited_old_clusters(old_cluster_to_paths_index.size(), false);

    for (uint32_t i = 0; i < old_cluster_to_paths_index.size(); ++i) {

        if (!visited_old_clusters.at(i)) {

            queue<uint32_t> search_queue;
            search_queue.push(i);

            cluster_to_paths_index.emplace_back(vector<uint32_t>());

            while (!search_queue.empty()) {

                auto cur_cluster = search_queue.front();
                
                if (!visited_old_clusters.at(cur_cluster)) {

                    visited_old_clusters.at(cur_cluster) = true;
                    cluster_to_paths_index.back().insert(cluster_to_paths_index.back().end(), old_cluster_to_paths_index.at(cur_cluster).begin(), old_cluster_to_paths_index.at(cur_cluster).end());

                    for (auto & next_cluster: connected_clusters.at(cur_cluster)) {

                        if (!visited_old_clusters.at(next_cluster)) {

                            search_queue.push(next_cluster);
                        }
                    }
                }

                search_queue.pop();
            }

            sort(cluster_to_paths_index.back().begin(), cluster_to_paths_index.back().end());

            for (auto & path: cluster_to_paths_index.back()) {

                path_to_cluster_index.at(path) = cluster_to_paths_index.size() - 1;
            } 
        }
    } 
}

