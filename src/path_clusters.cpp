#include <assert.h>
#include <queue>
#include <algorithm>

#include "path_clusters.hpp"
#include "utils.hpp"


static const uint32_t mutex_freq = 100;

PathClusters::PathClusters(const uint32_t num_threads_in, const PathsIndex & paths_index, const spp::sparse_hash_map<vector<AlignmentPath>, uint32_t> & align_paths_index) : num_threads(num_threads_in), num_paths(paths_index.numberOfPaths()) {

    vector<spp::sparse_hash_set<uint32_t> > connected_paths(num_paths, spp::sparse_hash_set<uint32_t>());
    vector<mutex> connected_paths_mutexes(ceil(num_paths / static_cast<double>(mutex_freq)));

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

                uint32_t anchor_path_id = numeric_limits<uint32_t>::max();

                for (size_t j = 0; j < align_paths_index_it->first.size() - 1; ++j) {

                    auto align_path_ids = paths_index.locatePathIds(align_paths_index_it->first.at(j).gbwt_search);
                    assert(!align_path_ids.empty());

                    if (anchor_path_id == numeric_limits<uint32_t>::max()) {

                        anchor_path_id = align_path_ids.front();
                    }

                    for (auto & path_id: align_path_ids) {

                        addConnection(&connected_paths, &connected_paths_mutexes, anchor_path_id, path_id);
                    }
                }

                advance(align_paths_index_it, num_threads);
                cur_align_paths_idx += num_threads;
            }
        }
    }

    createPathClusters(&path_to_cluster_index, &cluster_to_paths_index, connected_paths);
}

void PathClusters::addNodeClustering(const PathsIndex & paths_index, const bool cluster_strands) {

    auto node_connected_clusters = findNodeConnectedClusters(paths_index, cluster_strands);
    mergeClusters(&path_to_cluster_index, &cluster_to_paths_index, node_connected_clusters);
}

void PathClusters::createNoteClusterIndex(const PathsIndex & paths_index, const bool cluster_strands) {

    auto node_connected_paths = findNodeConnectedPaths(paths_index, cluster_strands);

    vector<vector<uint32_t> > note_cluster_to_paths_index;
    createPathClusters(&path_to_note_cluster_index, &note_cluster_to_paths_index, node_connected_paths);
}

void PathClusters::addConnection(vector<spp::sparse_hash_set<uint32_t> > * connections, vector<mutex> * connections_mutexes, const uint32_t anchor_id, const uint32_t id) const {

    if (anchor_id != id) {

        auto anchor_mutex_idx = floor(anchor_id / static_cast<double>(mutex_freq));
        auto mutex_idx = floor(id / static_cast<double>(mutex_freq));

        if (anchor_mutex_idx != mutex_idx) {

            lock(connections_mutexes->at(anchor_mutex_idx), connections_mutexes->at(mutex_idx));

        } else {

            connections_mutexes->at(anchor_mutex_idx).lock();
        }
                                
        connections->at(anchor_id).emplace(id);
        connections->at(id).emplace(anchor_id);

        connections_mutexes->at(anchor_mutex_idx).unlock();

        if (anchor_mutex_idx != mutex_idx) {

            connections_mutexes->at(mutex_idx).unlock();
        }
    }
}

vector<spp::sparse_hash_set<uint32_t> > PathClusters::findNodeConnectedPaths(const PathsIndex & paths_index, const bool cluster_strands) const {

    vector<spp::sparse_hash_set<uint32_t> > connected_paths(num_paths, spp::sparse_hash_set<uint32_t>());
    vector<mutex> connected_paths_mutexes(ceil(num_paths / static_cast<double>(mutex_freq)));

    #pragma omp parallel num_threads(num_threads)
    {
        #pragma omp for schedule(static)
        for (size_t i = 1; i <= paths_index.numberOfNodes(); ++i) {

            auto node_path_id_sets = getNodePathSets(i, paths_index);
            assert(node_path_id_sets.size() <= 2);

            uint32_t anchor_path_id = numeric_limits<uint32_t>::max();

            for (auto & node_path_ids: node_path_id_sets) {

                if (!node_path_ids.empty()) {

                    if (anchor_path_id == numeric_limits<uint32_t>::max() || !cluster_strands) {

                        anchor_path_id = node_path_ids.front();
                    }

                    for (auto & path_id: node_path_ids) {

                        addConnection(&connected_paths, &connected_paths_mutexes, anchor_path_id, path_id);
                    }
                }
            }
        }
    }

    return connected_paths;
}

vector<spp::sparse_hash_set<uint32_t> > PathClusters::findNodeConnectedClusters(const PathsIndex & paths_index, const bool cluster_strands) const {

    vector<spp::sparse_hash_set<uint32_t> > connected_clusters(cluster_to_paths_index.size(), spp::sparse_hash_set<uint32_t>());
    vector<mutex> connected_clusters_mutexes(ceil(cluster_to_paths_index.size() / static_cast<double>(mutex_freq)));

    #pragma omp parallel num_threads(num_threads)
    {
        #pragma omp for schedule(static)
        for (size_t i = 1; i <= paths_index.numberOfNodes(); ++i) {

            auto node_path_id_sets = getNodePathSets(i, paths_index);
            assert(node_path_id_sets.size() <= 2);

            uint32_t anchor_path_id = numeric_limits<uint32_t>::max();

            for (auto & node_path_ids: node_path_id_sets) {

                if (!node_path_ids.empty()) {

                    if (anchor_path_id == numeric_limits<uint32_t>::max() || !cluster_strands) {

                        anchor_path_id = node_path_ids.front();
                    }

                    for (auto & path_id: node_path_ids) {

                        addConnection(&connected_clusters, &connected_clusters_mutexes, path_to_cluster_index.at(anchor_path_id), path_to_cluster_index.at(path_id));
                    }
                }
            }
        }
    }

    return connected_clusters;
}

vector<vector<gbwt::size_type> > PathClusters::getNodePathSets(const size_t node_idx, const PathsIndex & paths_index) const {

    vector<vector<gbwt::size_type> > node_path_id_sets;

    pair<gbwt::SearchState, gbwt::size_type> gbwt_search;
    paths_index.find(&gbwt_search, gbwt::Node::encode(node_idx, false));

    if (!gbwt_search.first.empty()) {

        node_path_id_sets.emplace_back(paths_index.locatePathIds(gbwt_search));
    }

    if (!paths_index.bidirectional()) {

        paths_index.find(&gbwt_search, gbwt::Node::encode(node_idx, true));

        if (!gbwt_search.first.empty()) {

            node_path_id_sets.emplace_back(paths_index.locatePathIds(gbwt_search));
        }
    }

    return node_path_id_sets;
}

void PathClusters::createPathClusters(vector<uint32_t> * path_to_cluster, vector<vector<uint32_t> > * cluster_to_paths, const vector<spp::sparse_hash_set<uint32_t> > & connected_paths) {

    assert(path_to_cluster->empty());
    assert(cluster_to_paths->empty());

    *path_to_cluster = vector<uint32_t>(num_paths, -1);

    for (uint32_t i = 0; i < num_paths; ++i) {

        if (path_to_cluster->at(i) == -1) {

            queue<uint32_t> search_queue;
            search_queue.push(i);

            cluster_to_paths->emplace_back(vector<uint32_t>());

            while (!search_queue.empty()) {

                auto cur_path = search_queue.front();

                bool is_first_visit = (path_to_cluster->at(cur_path) == -1);
                assert(is_first_visit || path_to_cluster->at(cur_path) == cluster_to_paths->size() - 1);

                path_to_cluster->at(cur_path) = cluster_to_paths->size() - 1;
                
                if (is_first_visit) {

                    cluster_to_paths->back().emplace_back(cur_path);

                    for (auto & next_path: connected_paths.at(cur_path)) {

                        if (path_to_cluster->at(next_path) == -1) {

                            search_queue.push(next_path);
                        }
                    }
                }

                search_queue.pop();
            }

            sort(cluster_to_paths->back().begin(), cluster_to_paths->back().end());
        }
    }
}

void PathClusters::mergeClusters(vector<uint32_t> * path_to_cluster, vector<vector<uint32_t> > * cluster_to_paths, const vector<spp::sparse_hash_set<uint32_t> > & connected_clusters) {

    auto old_cluster_to_paths = *cluster_to_paths;
    cluster_to_paths->clear();

    vector<bool> visited_old_clusters(old_cluster_to_paths.size(), false);

    for (uint32_t i = 0; i < old_cluster_to_paths.size(); ++i) {

        if (!visited_old_clusters.at(i)) {

            queue<uint32_t> search_queue;
            search_queue.push(i);

            cluster_to_paths->emplace_back(vector<uint32_t>());

            while (!search_queue.empty()) {

                auto cur_cluster = search_queue.front();
                
                if (!visited_old_clusters.at(cur_cluster)) {

                    visited_old_clusters.at(cur_cluster) = true;
                    cluster_to_paths->back().insert(cluster_to_paths->back().end(), old_cluster_to_paths.at(cur_cluster).begin(), old_cluster_to_paths.at(cur_cluster).end());

                    for (auto & next_cluster: connected_clusters.at(cur_cluster)) {

                        if (!visited_old_clusters.at(next_cluster)) {

                            search_queue.push(next_cluster);
                        }
                    }
                }

                search_queue.pop();
            }

            sort(cluster_to_paths->back().begin(), cluster_to_paths->back().end());

            for (auto & path: cluster_to_paths->back()) {

                path_to_cluster->at(path) = cluster_to_paths->size() - 1;
            } 
        }
    } 
}

