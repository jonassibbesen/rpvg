
#include <assert.h>
#include <queue>
#include <algorithm>

#include "gbwt/gbwt.h"

#include "path_clusters.hpp"
#include "utils.hpp"


void PathClusters::findPathClusters(spp::sparse_hash_map<uint32_t, spp::sparse_hash_set<uint32_t> > * connected_paths, const PathsIndex & paths_index, const bool use_path_node_clustering) {

    if (use_path_node_clustering) {

        addPathNodeClusters(connected_paths, paths_index);
    }

    createPathClusters(*connected_paths, paths_index.index().metadata.paths());
}

void PathClusters::findPathClusters(vector<spp::sparse_hash_map<uint32_t, spp::sparse_hash_set<uint32_t> > > * connected_paths, const PathsIndex & paths_index, const bool use_path_node_clustering) {

    for (size_t i = 1; i < connected_paths->size(); ++i) {

        for (auto & path_component: connected_paths->at(i)) {

            auto connected_paths_it = connected_paths->front().emplace(path_component.first, spp::sparse_hash_set<uint32_t>());

            for (auto & path_id: path_component.second) {

                connected_paths_it.first->second.emplace(path_id);
            }
        }
    }

    findPathClusters(&(connected_paths->front()), paths_index, use_path_node_clustering);
}

void PathClusters::addPathNodeClusters(spp::sparse_hash_map<uint32_t, spp::sparse_hash_set<uint32_t> > * connected_paths, const PathsIndex & paths_index) {

    spp::sparse_hash_map<uint32_t, spp::sparse_hash_set<uint32_t> > connected_nodes;
    uint32_t max_node_id = 0;

    for (size_t i = 0; i < paths_index.index().sequences(); ++i) {

        auto gbwt_path = paths_index.index().extract(i);

        assert(!gbwt_path.empty());
        auto anchor_node_id = gbwt::Node::id(gbwt_path.front());

        for (auto & gbwt_node: gbwt_path) {

            uint32_t node_id = gbwt::Node::id(gbwt_node);
            max_node_id = max(max_node_id, node_id);

            if (anchor_node_id != node_id) {

                auto connected_nodes_it = connected_nodes.emplace(anchor_node_id, spp::sparse_hash_set<uint32_t>());
                connected_nodes_it.first->second.emplace(node_id);

                connected_nodes_it = connected_nodes.emplace(node_id, spp::sparse_hash_set<uint32_t>());
                connected_nodes_it.first->second.emplace(anchor_node_id);
            }
        }
    }

    PathClusters node_clusters;
    node_clusters.createPathClusters(connected_nodes, max_node_id + 1);
    
    unordered_map<uint32_t, vector<uint32_t> > path_clusters;

    for (size_t i = 0; i < paths_index.index().sequences(); ++i) {

        auto anchor_node_id = gbwt::Node::id(paths_index.index().extract(i).front());

        auto path_id = i;

        if (paths_index.index().bidirectional()) {

            path_id = gbwt::Path::id(i);
        }

        auto path_clusters_it = path_clusters.emplace(node_clusters.path_to_cluster_index.at(anchor_node_id), vector<uint32_t>());
        path_clusters_it.first->second.emplace_back(path_id);
    }

    for (auto & path_cluster: path_clusters) {

        assert(!path_cluster.second.empty());
        auto anchor_path_id = path_cluster.second.front();

        for (auto & path_id: path_cluster.second) {

            if (anchor_path_id != path_id) {

                auto connected_paths_it = connected_paths->emplace(anchor_path_id, spp::sparse_hash_set<uint32_t>());
                connected_paths_it.first->second.emplace(path_id);

                connected_paths_it = connected_paths->emplace(path_id, spp::sparse_hash_set<uint32_t>());
                connected_paths_it.first->second.emplace(anchor_path_id);
            }
        }
    }
}


void PathClusters::createPathClusters(const spp::sparse_hash_map<uint32_t, spp::sparse_hash_set<uint32_t> > & connected_paths, const uint32_t num_paths) {

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

            sort(cluster_to_paths_index.back().begin(), cluster_to_paths_index.back().end());
        }
    }
}
