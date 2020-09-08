#include <assert.h>
#include <queue>
#include <algorithm>

#include "path_clusters.hpp"
#include "utils.hpp"


PathClusters::PathClusters(const PathsIndex & paths_index_in, const uint32_t num_threads_in) : paths_index(paths_index_in), num_threads(num_threads_in) {

    connected_paths_t connected_paths;

    #pragma omp parallel num_threads(num_threads)
    {
        connected_paths_t thread_connected_paths;
        spp::sparse_hash_map<uint32_t, uint32_t> thread_node_to_path_index;

        #pragma omp for schedule(static)
        for (size_t i = 1; i <= paths_index.numberOfNodes(); ++i) {

            auto gbwt_search = paths_index.index().find(gbwt::Node::encode(i, false));
            vector<gbwt::size_type> node_path_ids;

            if (!gbwt_search.empty()) {

                node_path_ids = paths_index.locatePathIds(gbwt_search);
            }

            if (!paths_index.index().bidirectional()) {

                gbwt_search = paths_index.index().find(gbwt::Node::encode(i, true));

                if (!gbwt_search.empty()) {

                    auto node_path_ids_rev = paths_index.locatePathIds(gbwt_search);
                    node_path_ids.insert(node_path_ids.end(), node_path_ids_rev.begin(), node_path_ids_rev.end());
                }
            }

            if (!node_path_ids.empty()) {

                auto anchor_path_id = node_path_ids.front();

                for (auto & path_id: node_path_ids) {

                    if (anchor_path_id != path_id) {

                        auto thread_connected_paths_it = thread_connected_paths.emplace(anchor_path_id, spp::sparse_hash_set<uint32_t>());
                                                
                        if (thread_connected_paths_it.first->second.emplace(path_id).second) {

                            thread_connected_paths_it = thread_connected_paths.emplace(path_id, spp::sparse_hash_set<uint32_t>());
                            thread_connected_paths_it.first->second.emplace(anchor_path_id);
                        }
                    }
                }

                thread_node_to_path_index.emplace(i, anchor_path_id);
            }
        }

        #pragma omp critical
        {

            cerr << thread_connected_paths.size() << " " << thread_node_to_path_index.size() << endl;

            for (auto & paths: thread_connected_paths) {

                auto connected_paths_it = connected_paths.emplace(paths.first, spp::sparse_hash_set<uint32_t>());
                connected_paths_it.first->second.insert(paths.second.begin(), paths.second.end());
            }

            for (auto & node_path: thread_node_to_path_index) {

                assert(node_to_path_index.emplace(node_path).second);
            }
        }
    }

    createPathClusters(connected_paths);
}

void PathClusters::addReadClusters(const spp::sparse_hash_map<vector<AlignmentPath>, uint32_t> & align_paths_index) {

    connected_paths_t multimap_connected_clusters;

    #pragma omp parallel num_threads(num_threads)
    {
        uint32_t num_multimap = 0;
        connected_paths_t threaded_multimap_connected_clusters;

        #pragma omp for schedule(static, 1)
        for (size_t i = 0; i < num_threads; ++i) {

            uint32_t cur_align_paths_index_pos = 0;

            for (auto & align_paths: align_paths_index) {

                assert(!align_paths.first.empty());

                if (align_paths.first.front().is_multimap && (cur_align_paths_index_pos % num_threads == i)) { 

                    num_multimap++;
                    uint32_t anchor_path_cluster = 0;

                    for (size_t j = 0; j < align_paths.first.size(); ++j) {

                        assert(align_paths.first.at(j).is_multimap);
                        auto align_path_ids = paths_index.locatePathIds(align_paths.first.at(j).search_state);

                        if (j == 0) {

                            anchor_path_cluster = path_to_cluster_index.at(align_path_ids.front());
                        }

                        for (auto & path_id: align_path_ids) {

                            auto path_cluster = path_to_cluster_index.at(path_id);

                            if (anchor_path_cluster != path_cluster) {
                                
                                auto threaded_multimap_connected_clusters_it = threaded_multimap_connected_clusters.emplace(anchor_path_cluster, spp::sparse_hash_set<uint32_t>());
                                
                                if (threaded_multimap_connected_clusters_it.first->second.emplace(path_cluster).second) {

                                    threaded_multimap_connected_clusters_it = threaded_multimap_connected_clusters.emplace(path_cluster, spp::sparse_hash_set<uint32_t>());
                                    threaded_multimap_connected_clusters_it.first->second.emplace(anchor_path_cluster);
                                }
                            }
                        }
                    }
                }

                ++cur_align_paths_index_pos;
            }
        }

        #pragma omp critical
        {
            cerr << num_multimap << " " << threaded_multimap_connected_clusters.size() << endl;

            for (auto & paths: threaded_multimap_connected_clusters) {

                auto multimap_connected_clusters_it = multimap_connected_clusters.emplace(paths.first, spp::sparse_hash_set<uint32_t>());
                multimap_connected_clusters_it.first->second.insert(paths.second.begin(), paths.second.end());
            }
        }
    }

    mergeClusters(multimap_connected_clusters);
}

// void PathClusters::addCallTraversalClusters() {

//     spp::sparse_hash_map<string, vector<uint32_t> > connected_path_names;

//     for (size_t i = 0; i < paths_index.index().sequences(); ++i) {

//         auto path_id = i;

//         if (paths_index.index().bidirectional()) {

//             path_id = gbwt::Path::id(i);
//         }

//         auto call_path_name_split = splitString(paths_index.pathName(path_id), '_');
//         assert(call_path_name_split.size() == 4);

//         stringstream call_path_name_ss;
//         call_path_name_ss << call_path_name_split.at(0);
//         call_path_name_ss << "_" << call_path_name_split.at(1);
//         call_path_name_ss << "_" << call_path_name_split.at(2);

//         auto connected_path_names_it = connected_path_names.emplace(call_path_name_ss.str(), vector<uint32_t>());
//         connected_path_names_it.first->second.emplace_back(path_id);
//     }

//     auto connected_paths = findConnectedPaths();

//     for (auto & path_names: connected_path_names) {

//         assert(!path_names.second.empty());
//         auto anchor_path_id = path_names.second.front();

//         for (auto & path_id: path_names.second) {

//             if (anchor_path_id != path_id) {

//                 auto connected_paths_it = connected_paths.emplace(anchor_path_id, spp::sparse_hash_set<uint32_t>());
                
//                 if (connected_paths_it.first->second.emplace(path_id).second) {

//                     connected_paths_it = connected_paths.emplace(path_id, spp::sparse_hash_set<uint32_t>());
//                     connected_paths_it.first->second.emplace(anchor_path_id);
//                 }
//             }
//         }
//     } 

//     createPathClusters(connected_paths);
// }

void PathClusters::createPathClusters(const connected_paths_t & connected_paths) {

    path_to_cluster_index.clear();
    cluster_to_paths_index.clear();

    path_to_cluster_index = vector<uint32_t>(paths_index.index().metadata.paths(), -1);

    for (uint32_t i = 0; i < paths_index.index().metadata.paths(); ++i) {

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

                        for (auto & next_path: connected_paths_it->second) {

                            if (path_to_cluster_index.at(next_path) == -1) {

                                search_queue.push(next_path);
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

void PathClusters::mergeClusters(const connected_paths_t & connected_clusters) {

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

                    auto connected_clusters_it = connected_clusters.find(cur_cluster);

                    if (connected_clusters_it != connected_clusters.end()) {

                        for (auto & next_cluster: connected_clusters_it->second) {

                            if (!visited_old_clusters.at(next_cluster)) {

                                search_queue.push(next_cluster);
                            }
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

