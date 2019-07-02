
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <chrono>
#include <limits>

#include <cxxopts.hpp>
#include <gbwt/gbwt.h>
#include <gssw.h>
#include <vg/io/vpkg.hpp>
#include <vg/io/stream.hpp>
#include <vg/io/basic_stream.hpp>

#include "io/register_libvg_io.hpp"

#include "utils.hpp"
#include "alignment_path.hpp"
#include "path_clusters.hpp"
#include "read_path_probs.hpp"

using namespace std;

const double frag_length_mean = 277;
const double frag_length_sd = 43;

// #define debug


vector<AlignmentPath> get_align_paths(const vg::Alignment & alignment, const gbwt::GBWT & paths_index) {
            
    AlignmentPath align_path;
    align_path.scores.first = alignment.score();
    align_path.mapqs.first = alignment.mapping_quality();

    align_path.extendPath(alignment.path(), 0, paths_index);

    return vector<AlignmentPath>(1, align_path);
}

vector<AlignmentPath> get_align_paths(const vg::MultipathAlignment & mp_alignment, const gbwt::GBWT & paths_index) {

    vector<AlignmentPath> align_paths;

    std::queue<pair<AlignmentPath, int32_t> > mp_align_path_queue;
    assert(mp_alignment.start_size() > 0); 

    for (auto & start_subpath_idx: mp_alignment.start()) {

        mp_align_path_queue.push(make_pair(AlignmentPath(), start_subpath_idx));
    }

    // Perform depth-first alignment path extension.
    while (!mp_align_path_queue.empty()) {

        auto & cur_mp_align_path = mp_align_path_queue.front();

        const vg::Subpath & subpath = mp_alignment.subpath(cur_mp_align_path.second);

        cur_mp_align_path.first.scores.first += subpath.score();
        cur_mp_align_path.first.extendPath(subpath.path(), 0, paths_index);

        if (subpath.next_size() > 0) {

            for (auto & next_subpath_idx: subpath.next()) {

                mp_align_path_queue.push(make_pair(cur_mp_align_path.first, next_subpath_idx));
            }

        } else {

            align_paths.emplace_back(cur_mp_align_path.first);
        }

        mp_align_path_queue.pop();
    }

    return align_paths;
}

template<typename AlignType>
vector<AlignmentPath> get_align_paths_with_ids(const AlignType & alignment, const gbwt::GBWT & paths_index) {

    auto align_paths = get_align_paths(alignment, paths_index);

    for (auto & align_path: align_paths) {

        align_path.path_ids = paths_index.locate(align_path.path);
    }

    return align_paths;
}

multimap<gbwt::node_type, pair<int32_t, int32_t> > get_gbwt_node_index(const vg::Alignment & alignment) {

    multimap<gbwt::node_type, pair<int32_t, int32_t> > gbwt_node_index;

    for (size_t i = 0; i < alignment.path().mapping_size(); ++i) {

        gbwt_node_index.emplace(mapping_to_gbwt(alignment.path().mapping(i)), make_pair(0, i));
    }

    return gbwt_node_index;
}

void find_paired_align_paths(vector<AlignmentPath> * paired_align_paths, const AlignmentPath & start_align_path, const vg::Alignment & end_alignment, const gbwt::GBWT & paths_index, const vector<uint32_t> & node_seq_lengths, const int32_t max_pair_distance) {

    assert(!start_align_path.path.empty());

    std::queue<AlignmentPath> paired_align_path_queue;
    paired_align_path_queue.push(start_align_path);

    assert(end_alignment.path().mapping_size() > 0);

    auto end_alignment_node_index = get_gbwt_node_index(end_alignment);

    // Perform depth-first path extension.
    while (!paired_align_path_queue.empty()) {

        auto & cur_paired_align_path = paired_align_path_queue.front();

        if (cur_paired_align_path.seq_length > max_pair_distance) {

            paired_align_path_queue.pop();
            continue;                
        }

        auto end_alignment_node_index_it = end_alignment_node_index.equal_range(cur_paired_align_path.path.node);

        // Stop current extension if end node is reached.
        if (end_alignment_node_index_it.first != end_alignment_node_index_it.second) {

            while (end_alignment_node_index_it.first != end_alignment_node_index_it.second) {

                auto start_mapping_idx = end_alignment_node_index_it.first->second.second + 1;

                for (size_t i = start_mapping_idx; i < end_alignment.path().mapping_size(); ++i) {

                    auto gbwt_node_id = mapping_to_gbwt(end_alignment.path().mapping(i));

                    cur_paired_align_path.path = paths_index.extend(cur_paired_align_path.path, gbwt_node_id);
                    cur_paired_align_path.node_length++;
                    cur_paired_align_path.seq_length += node_seq_lengths.at(gbwt::Node::id(gbwt_node_id));
                }

                cur_paired_align_path.seq_length -= (node_seq_lengths.at(end_alignment.path().mapping().rbegin()->position().node_id()) - mapping_to_length(*end_alignment.path().mapping().rbegin()));

                if (!cur_paired_align_path.path.empty() && cur_paired_align_path.seq_length <= max_pair_distance) {

                    paired_align_paths->emplace_back(cur_paired_align_path);

                    paired_align_paths->back().path_ids = paths_index.locate(paired_align_paths->back().path);                                
                    paired_align_paths->back().scores.second = end_alignment.score();
                    paired_align_paths->back().mapqs.second = end_alignment.mapping_quality();           
                }

                end_alignment_node_index_it.first++;
            }

            paired_align_path_queue.pop();
            continue;
        }
        
        auto out_edges = paths_index.edges(cur_paired_align_path.path.node);

        // End current extension if no outgoing edges exist.
        if (out_edges.empty()) {

            paired_align_path_queue.pop();
            continue;
        }

        auto out_edges_it = out_edges.begin(); 

        while (out_edges_it != out_edges.end()) {

            if (out_edges_it->first != gbwt::ENDMARKER) {

                auto extended_path = paths_index.extend(cur_paired_align_path.path, out_edges_it->first);

                // Add new extension to queue if not empty (path found).
                if (!extended_path.empty()) { 

                    AlignmentPath new_paired_align_path = cur_paired_align_path;
                    new_paired_align_path.path = extended_path;
                    new_paired_align_path.node_length++;
                    new_paired_align_path.seq_length += node_seq_lengths.at(gbwt::Node::id(out_edges_it->first));

                    paired_align_path_queue.push(new_paired_align_path);
                }
            }

            ++out_edges_it;
        }

        paired_align_path_queue.pop();
    }
}

vector<AlignmentPath> get_paired_align_paths(const vg::Alignment & alignment_1, const vg::Alignment & alignment_2, const gbwt::GBWT & paths_index, const vector<uint32_t> & node_seq_lengths, const int32_t max_pair_distance) {

    vector<AlignmentPath> paired_align_paths;

    function<size_t(const int64_t)> node_seq_length_func = [&node_seq_lengths](const int64_t node_id) { return node_seq_lengths.at(node_id); };

    auto align_path_1 = get_align_paths(alignment_1, paths_index);
    assert(align_path_1.size() == 1);

    if (!align_path_1.front().path.empty()) {

        vg::Alignment alignment_2_rc = alignment_2; 
        *alignment_2_rc.mutable_path() = lazy_reverse_complement_path(alignment_2.path(), node_seq_length_func);
        find_paired_align_paths(&paired_align_paths, align_path_1.front(), alignment_2_rc, paths_index, node_seq_lengths, max_pair_distance);
    }

    auto align_path_2 = get_align_paths(alignment_2, paths_index);
    assert(align_path_2.size() == 1);

    if (!align_path_2.front().path.empty()) {

        vg::Alignment alignment_1_rc = alignment_1; 
        *alignment_1_rc.mutable_path() = lazy_reverse_complement_path(alignment_1.path(), node_seq_length_func);
        find_paired_align_paths(&paired_align_paths, align_path_2.front(), alignment_1_rc, paths_index, node_seq_lengths, max_pair_distance);
    }

    return paired_align_paths;
}

int main(int argc, char* argv[]) {

    cxxopts::Options options("vgprob", "calculate read path probabilities");

    options.add_options()
      ("g,graph", "vg graph file name (required)", cxxopts::value<string>())
      ("p,paths", "gbwt index file name (required)", cxxopts::value<string>())
      ("a,alignments", "gam alignments file name (required)", cxxopts::value<string>())
      ("o,output", "output prefix (required)", cxxopts::value<string>())
      // ("m,multipath", "alignment input is multipath (gamp)", cxxopts::value<bool>())
      ("t,threads", "number of compute threads", cxxopts::value<int32_t>()->default_value("1"))
      ("h,help", "print help", cxxopts::value<bool>())
      ;

    if (argc == 1) {

        cout << options.help() << endl;
        return 1;
    }

    auto option_results = options.parse(argc, argv);

    if (option_results.count("help")) {

        cout << options.help() << endl;
        return 1;
    }

    assert(option_results.count("graph") == 1);
    assert(option_results.count("paths") == 1);
    assert(option_results.count("alignments") == 1);
    assert(option_results.count("output") == 1);
    assert(option_results.count("threads") <= 1);

    const int32_t num_threads = option_results["threads"].as<int32_t>();

    assert(num_threads > 0);
    omp_set_num_threads(num_threads);

    double time1 = gbwt::readTimer();

    vector<uint32_t> node_seq_lengths;

    {
        vg::Graph graph = vg::io::inputStream(option_results["graph"].as<string>());
        node_seq_lengths = vector<uint32_t>(graph.node_size() + 1, 0);

        for (auto & node: graph.node()) {

            if (node.id() >= node_seq_lengths.size()) {

                node_seq_lengths.resize(node.id() + 1, 0);
            }

            assert(node_seq_lengths.at(node.id()) == 0);
            node_seq_lengths.at(node.id()) = node.sequence().size();
        }
    }

    double time2 = gbwt::readTimer();
    cout << "Load graph and find node seq_lengths " << time2 - time1 << " seconds, " << gbwt::inGigabytes(gbwt::memoryUsage()) << " GB" << endl;

    assert(vg::io::register_libvg_io());

    unique_ptr<gbwt::GBWT> paths_index = vg::io::VPKG::load_one<gbwt::GBWT>(option_results["paths"].as<string>());

    double time3 = gbwt::readTimer();
    cout << "Load GBWT " << time3 - time2 << " seconds, " << gbwt::inGigabytes(gbwt::memoryUsage()) << " GB" << endl;

    ifstream alignment_istream(option_results["alignments"].as<string>());
    assert(alignment_istream.is_open());

    ofstream probs_ostream(option_results["output"].as<string>() + ".txt");
    assert(probs_ostream.is_open());

    const auto num_paths = paths_index->metadata.haplotype_count;

    vector<unordered_map<int32_t, unordered_set<int32_t> > > connected_paths_threads(num_threads);
    vector<vector<vector<AlignmentPath> > > paired_align_paths_threads(num_threads);

    if (option_results.count("multipath")) {

        vg::io::for_each_interleaved_pair_parallel<vg::MultipathAlignment>(alignment_istream, [&](vg::MultipathAlignment & mp_alignment_1, vg::MultipathAlignment & mp_alignment_2) {

            if (mp_alignment_1.subpath_size() > 0 && mp_alignment_2.subpath_size() > 0) {

    #ifdef debug
                cout << endl;
                cout << pb2json(mp_alignment_1) << endl;
                cout << pb2json(mp_alignment_2) << endl;


                vg::MultipathAlignment mp_alignment_1_rc = mp_alignment_1; 
                mp_alignment_1_rc.clear_start();
                *mp_alignment_1_rc.mutable_subpath() = lazy_reverse_complement_subpaths(mp_alignment_1_rc.subpath(), mp_alignment_1_rc.mutable_start(), [&node_seq_lengths](const int64_t node_id) { return node_seq_lengths.at(node_id); });

                vg::MultipathAlignment mp_alignment_2_rc = mp_alignment_2; 
                mp_alignment_2_rc.clear_start();
                *mp_alignment_2_rc.mutable_subpath() = lazy_reverse_complement_subpaths(mp_alignment_2_rc.subpath(), mp_alignment_2_rc.mutable_start(), [&node_seq_lengths](const int64_t node_id) { return node_seq_lengths.at(node_id); });

                cout << pb2json(mp_alignment_1_rc) << endl;
                cout << pb2json(mp_alignment_2_rc) << endl;

                cout << get_align_paths_with_ids<vg::MultipathAlignment>(mp_alignment_1, *paths_index) << endl;
                cout << get_align_paths_with_ids<vg::MultipathAlignment>(mp_alignment_1_rc, *paths_index) << endl;
                cout << get_align_paths_with_ids<vg::MultipathAlignment>(mp_alignment_2, *paths_index) << endl;
                cout << get_align_paths_with_ids<vg::MultipathAlignment>(mp_alignment_2_rc, *paths_index) << endl;       
    #endif 

            }
        });

    } else {

        vg::io::for_each_interleaved_pair_parallel<vg::Alignment>(alignment_istream, [&](vg::Alignment & alignment_1, vg::Alignment & alignment_2) {

            if (alignment_1.has_path() && alignment_2.has_path()) {

                auto paired_align_paths = get_paired_align_paths(alignment_1, alignment_2, *paths_index, node_seq_lengths, frag_length_mean + 10 * frag_length_sd);

    #ifdef debug
                cout << endl;
                cout << pb2json(alignment_1) << endl;
                cout << pb2json(alignment_2) << endl;

                vg::Alignment alignment_1_rc = alignment_1; 
                *alignment_1_rc.mutable_path() = lazy_reverse_complement_path(alignment_1.path(), [&node_seq_lengths](const int64_t node_id) { return node_seq_lengths.at(node_id); });

                vg::Alignment alignment_2_rc = alignment_2; 
                *alignment_2_rc.mutable_path() = lazy_reverse_complement_path(alignment_2.path(), [&node_seq_lengths](const int64_t node_id) { return node_seq_lengths.at(node_id); });

                cout << pb2json(alignment_1_rc) << endl;
                cout << pb2json(alignment_2_rc) << endl;

                cout << get_align_paths_with_ids<vg::Alignment>(alignment_1, *paths_index) << endl;
                cout << get_align_paths_with_ids<vg::Alignment>(alignment_1_rc, *paths_index) << endl;
                cout << get_align_paths_with_ids<vg::Alignment>(alignment_2, *paths_index) << endl;
                cout << get_align_paths_with_ids<vg::Alignment>(alignment_2_rc, *paths_index) << endl;
          
                cout << paired_align_paths << endl;
    #endif 

                if (!paired_align_paths.empty()) {

                    auto anchor_path_id = paired_align_paths.front().path_ids.front();

                    for (auto & align_path: paired_align_paths) {

                        for (auto & path_id: align_path.path_ids) {

                            if (anchor_path_id != path_id) {

                                connected_paths_threads.at(omp_get_thread_num())[anchor_path_id].emplace(path_id);
                                connected_paths_threads.at(omp_get_thread_num())[path_id].emplace(anchor_path_id);
                            }
                        }
                    }

                    paired_align_paths_threads.at(omp_get_thread_num()).emplace_back(paired_align_paths);
                }
            }
        });
    }

    double time5 = gbwt::readTimer();
    cout << "Found paired alignment paths " << time5 - time3 << " seconds, " << gbwt::inGigabytes(gbwt::memoryUsage()) << " GB" << endl;

    for (size_t i = 1; i < connected_paths_threads.size(); ++i) {

        for (auto & connected_path_clusters: connected_paths_threads.at(i)) {

            auto connected_paths_threads_it = connected_paths_threads.front().emplace(connected_path_clusters.first, unordered_set<int32_t>());
            for (auto & paths: connected_path_clusters.second) {

                connected_paths_threads_it.first->second.emplace(paths);
            }
        }
    }

    double time6 = gbwt::readTimer();
    cout << "Merged connected path threads " << time6 - time5 << " seconds, " << gbwt::inGigabytes(gbwt::memoryUsage()) << " GB" << endl;

    auto path_clusters = PathClusters(connected_paths_threads.front(), num_paths);

    double time62 = gbwt::readTimer();
    cout << "Found path clusters " << time62 - time6 << " seconds, " << gbwt::inGigabytes(gbwt::memoryUsage()) << " GB" << endl;
 
    vector<vector<vector<AlignmentPath> > > clustered_paired_align_paths(path_clusters.cluster_to_path_index.size());

    for (auto & paired_align_paths: paired_align_paths_threads) {

        for (auto & align_path: paired_align_paths) {

            clustered_paired_align_paths.at(path_clusters.path_to_cluster_index.at(align_path.front().path_ids.front())).push_back(move(align_path));
        }
    }

    double time7 = gbwt::readTimer();
    cout << "Clustered paired alignment paths " << time7 - time6 << " seconds, " << gbwt::inGigabytes(gbwt::memoryUsage()) << " GB" << endl;

    vector<vector<ReadPathProbs> > clustered_paired_align_path_probs(path_clusters.cluster_to_path_index.size());

    #pragma omp parallel
    { 
    
        #pragma omp for
        for (size_t i = 0; i < clustered_paired_align_paths.size(); ++i) {

            clustered_paired_align_path_probs.at(i).reserve(clustered_paired_align_paths.at(i).size());

            unordered_map<uint32_t, uint32_t> clustered_path_index;

            for (auto & path_id: path_clusters.cluster_to_path_index.at(i)) {

                assert(clustered_path_index.emplace(path_id, clustered_path_index.size()).second);
            }

            for (auto & align_paths: clustered_paired_align_paths.at(i)) {

                clustered_paired_align_path_probs.at(i).emplace_back(clustered_path_index.size());
                clustered_paired_align_path_probs.at(i).back().calcReadPathProbs(align_paths, clustered_path_index, frag_length_mean, frag_length_sd);
            }

            sort(clustered_paired_align_path_probs.at(i).begin(), clustered_paired_align_path_probs.at(i).end());
        }
    }

    double time8 = gbwt::readTimer();
    cout << "Calculated and sorted paired alignment path probabilites " << time8 - time7 << " seconds, " << gbwt::inGigabytes(gbwt::memoryUsage()) << " GB" << endl;

    for (size_t i = 0; i < clustered_paired_align_path_probs.size(); ++i) {

        const vector<ReadPathProbs> & clustered_probs = clustered_paired_align_path_probs.at(i);

        if (clustered_probs.empty()) {

            continue;
        }

        probs_ostream << "#\nx Noise";
        for (auto & path_id: path_clusters.cluster_to_path_index.at(i)) {

            probs_ostream << " " << getPathName(*paths_index, path_id);
        }
        probs_ostream << endl;

        int32_t num_paired_paths = 1;
        int32_t prev_unique_probs_idx = 0;

        for (size_t j = 1; j < clustered_probs.size(); ++j) {

            if (clustered_probs.at(prev_unique_probs_idx) == clustered_probs.at(j)) {

                num_paired_paths++;
            
            } else {

                probs_ostream << num_paired_paths << " " << clustered_probs.at(prev_unique_probs_idx);
                num_paired_paths = 1;
                prev_unique_probs_idx = j;
            }
        }

        probs_ostream << num_paired_paths << " " << clustered_probs.at(prev_unique_probs_idx);
    }

    alignment_istream.close();
    probs_ostream.close();
 
    double time9 = gbwt::readTimer();
    cout << "Collapsed and wrote probabilites " << time9 - time8 << " seconds, " << gbwt::inGigabytes(gbwt::memoryUsage()) << " GB" << endl;

	return 0;
}