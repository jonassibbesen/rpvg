
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <chrono>
#include <limits>
#include <algorithm>

#include <cxxopts.hpp>
#include <gbwt/gbwt.h>
#include <vg/io/vpkg.hpp>
#include <vg/io/stream.hpp>
#include <vg/io/basic_stream.hpp>

#include "io/register_libvg_io.hpp"

#include "utils.hpp"
#include "alignment_path.hpp"
#include "alignment_path_finder.hpp"
#include "path_clusters.hpp"
#include "read_path_probs.hpp"

using namespace std;

const double frag_length_mean = 277;
const double frag_length_sd = 43;


void addPairedAlignmentPathsThreaded(vector<unordered_map<int32_t, unordered_set<int32_t> > > * connected_paths_threads, vector<vector<vector<AlignmentPath> > > * paired_align_paths_threads, const vector<AlignmentPath> & paired_align_paths, const int32_t thread_num) {

    if (!paired_align_paths.empty()) {

        auto anchor_path_id = paired_align_paths.front().path_ids.front();

        for (auto & align_path: paired_align_paths) {

            for (auto & path_id: align_path.path_ids) {

                if (anchor_path_id != path_id) {

                    connected_paths_threads->at(thread_num)[anchor_path_id].emplace(path_id);
                    connected_paths_threads->at(thread_num)[path_id].emplace(anchor_path_id);
                }
            }
        }

        paired_align_paths_threads->at(thread_num).emplace_back(paired_align_paths);
    }    
}


int main(int argc, char* argv[]) {

    cxxopts::Options options("vgprob", "calculate read path probabilities");

    options.add_options()
      ("g,graph", "vg graph file name (required)", cxxopts::value<string>())
      ("p,paths", "gbwt index file name (required)", cxxopts::value<string>())
      ("a,alignments", "gam alignments file name (required)", cxxopts::value<string>())
      ("o,output", "output prefix (required)", cxxopts::value<string>())
      ("m,multipath", "alignment input is multipath (gamp)", cxxopts::value<bool>())
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

    vg::Graph graph = vg::io::inputStream(option_results["graph"].as<string>());

    assert(vg::io::register_libvg_io());
    unique_ptr<gbwt::GBWT> paths_index = vg::io::VPKG::load_one<gbwt::GBWT>(option_results["paths"].as<string>());

    double time2 = gbwt::readTimer();
    cout << "Load graph and GBWT " << time2 - time1 << " seconds, " << gbwt::inGigabytes(gbwt::memoryUsage()) << " GB" << endl;

    ifstream alignment_istream(option_results["alignments"].as<string>());
    assert(alignment_istream.is_open());

    ofstream probs_ostream(option_results["output"].as<string>() + ".txt");
    assert(probs_ostream.is_open());

    const auto num_paths = paths_index->metadata.haplotype_count;

    vector<unordered_map<int32_t, unordered_set<int32_t> > > connected_paths_threads(num_threads);
    vector<vector<vector<AlignmentPath> > > paired_align_paths_threads(num_threads);

    if (option_results.count("multipath")) {

        AlignmentPathFinder<vg::MultipathAlignment> alignment_path_finder(graph, *paths_index);

        vg::io::for_each_interleaved_pair_parallel<vg::MultipathAlignment>(alignment_istream, [&](vg::MultipathAlignment & alignment_1, vg::MultipathAlignment & alignment_2) {

            if (alignment_1.subpath_size() > 0 && alignment_2.subpath_size() > 0) {

                auto paired_align_paths = alignment_path_finder.findPairedAlignmentPathsIds(alignment_1, alignment_2, frag_length_mean + 10 * frag_length_sd);
                addPairedAlignmentPathsThreaded(&connected_paths_threads, &paired_align_paths_threads, paired_align_paths, omp_get_thread_num());
            }
        });

    } else {

        AlignmentPathFinder<vg::Alignment> alignment_path_finder(graph, *paths_index);

        vg::io::for_each_interleaved_pair_parallel<vg::Alignment>(alignment_istream, [&](vg::Alignment & alignment_1, vg::Alignment & alignment_2) {

            if (alignment_1.has_path() && alignment_2.has_path()) {

                auto paired_align_paths = alignment_path_finder.findPairedAlignmentPathsIds(alignment_1, alignment_2, frag_length_mean + 10 * frag_length_sd);
                addPairedAlignmentPathsThreaded(&connected_paths_threads, &paired_align_paths_threads, paired_align_paths, omp_get_thread_num());
            }
        });
    }

    double time5 = gbwt::readTimer();
    cout << "Found paired alignment paths " << time5 - time2 << " seconds, " << gbwt::inGigabytes(gbwt::memoryUsage()) << " GB" << endl;

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

                probs_ostream << num_paired_paths << " " << clustered_probs.at(prev_unique_probs_idx) << endl;
                num_paired_paths = 1;
                prev_unique_probs_idx = j;
            }
        }

        probs_ostream << num_paired_paths << " " << clustered_probs.at(prev_unique_probs_idx) << endl;
    }

    alignment_istream.close();
    probs_ostream.close();
 
    double time9 = gbwt::readTimer();
    cout << "Collapsed and wrote probabilites " << time9 - time8 << " seconds, " << gbwt::inGigabytes(gbwt::memoryUsage()) << " GB" << endl;

	return 0;
}