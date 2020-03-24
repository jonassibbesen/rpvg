
#include <assert.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <chrono>
#include <limits>
#include <algorithm>
#include <iomanip>

#include "cxxopts.hpp"
#include "gbwt/gbwt.h"
#include "vg/io/vpkg.hpp"
#include "vg/io/stream.hpp"
#include "vg/io/basic_stream.hpp"
#include "io/register_libvg_io.hpp"
#include "handlegraph/handle_graph.hpp"
#include "gssw.h"
#include "utils.hpp"
#include "fragment_length_dist.hpp"
#include "paths_index.hpp"
#include "alignment_path.hpp"
#include "alignment_path_finder.hpp"
#include "path_clusters.hpp"
#include "read_path_probs.hpp"


void addAlignmentPathsThreaded(vector<unordered_map<int32_t, unordered_set<int32_t> > > * connected_paths_threads, vector<vector<vector<AlignmentPath> > > * align_paths_threads, const vector<AlignmentPath> & align_paths, const int32_t thread_idx) {

    if (!align_paths.empty()) {

        auto anchor_path_id = align_paths.front().ids.front();

        for (auto & align_path: align_paths) {

            for (auto & path_id: align_path.ids) {

                if (anchor_path_id != path_id) {

                    connected_paths_threads->at(thread_idx)[anchor_path_id].emplace(path_id);
                    connected_paths_threads->at(thread_idx)[path_id].emplace(anchor_path_id);
                }
            }
        }

        align_paths_threads->at(thread_idx).emplace_back(align_paths);
    }    
}

int main(int argc, char* argv[]) {

    cxxopts::Options options("rpvg", "calculate expression from read aligment path probabilities in variation graphs");

    options.add_options()
      ("g,graph", "xg graph file name (required)", cxxopts::value<string>())
      ("p,paths", "gbwt index file name (required)", cxxopts::value<string>())
      ("a,alignments", "gam(p) alignment file name (required)", cxxopts::value<string>())
      ("s,single-end", "alignment input is single-end reads", cxxopts::value<bool>())
      ("l,long-reads", "alignment input is non-fragmented long reads (single-end only)", cxxopts::value<bool>())
      ("u,multipath", "alignment input is multipath gamp format (default: gam)", cxxopts::value<bool>())
      ("o,output", "output file prefix", cxxopts::value<string>()->default_value("stdout"))
      ("m,frag-mean", "mean for fragment length distribution", cxxopts::value<double>())
      ("d,frag-sd", "standard deviation for fragment length distribution", cxxopts::value<double>())
      ("t,threads", "number of compute threads", cxxopts::value<int32_t>()->default_value("1"))
      ("h,help", "print help", cxxopts::value<bool>())
      ;

    if (argc == 1) {

        cerr << options.help() << endl;
        return 1;
    }

    auto option_results = options.parse(argc, argv);

    if (option_results.count("help")) {

        cerr << options.help() << endl;
        return 1;
    }

    assert(option_results.count("graph") == 1);
    assert(option_results.count("paths") == 1);
    assert(option_results.count("alignments") == 1);

    bool is_single_end = option_results.count("single-end");
    bool is_long_reads = option_results.count("long-reads");
    bool is_multipath = option_results.count("multipath");

    if (is_long_reads) {

        is_single_end = true;
    }

    if (option_results.count("frag-mean") != option_results.count("frag-sd")) {

        cerr << "ERROR: Both --frag-mean and --frag-sd needs to be given as input. Alternative, no values can be given for paired-end, non-long read alignments and the parameter estimated during mapping will be used instead (contained in the alignment file)." << endl;
        return 1;
    }

    FragmentLengthDist fragment_length_dist; 

    if (!option_results.count("frag-mean") && !option_results.count("frag-sd")) {

        if (is_single_end && !is_long_reads) {

            cerr << "ERROR: Both --frag-mean and --frag-sd needs to be given as input when using single-end, non-long read alignments." << endl;
            return 1;            
        }

        cerr << "WARNING: Fragment length distribution parameters not given as input. Parameters will be based on the first alignment that contains such information. Alternatively use --frag-mean and --frag-sd." << endl;

        ifstream frag_alignments_istream(option_results["alignments"].as<string>());
        assert(frag_alignments_istream.is_open());

        fragment_length_dist = FragmentLengthDist(&frag_alignments_istream, is_multipath);

        frag_alignments_istream.close();

        if (!fragment_length_dist.isValid()) {

            cerr << "ERROR: No fragment length distribution parameters found in alignments. Use --frag-mean and --frag-sd instead." << endl;
            return 1;
        
        } else {

            cerr << "Fragment length distribution parameters found in alignment (mean: " << fragment_length_dist.mean() << ", standard deviation: " << fragment_length_dist.sd() << ")" << endl;
        }      

    } else {

        fragment_length_dist = FragmentLengthDist(option_results["frag-mean"].as<double>(), option_results["frag-sd"].as<double>());
        cerr << "Fragment length distribution parameters given as input (mean: " << fragment_length_dist.mean() << ", standard deviation: " << fragment_length_dist.sd() << ")" << endl;
    }

    assert(fragment_length_dist.isValid());

    const int32_t num_threads = option_results["threads"].as<int32_t>();

    assert(num_threads > 0);
    omp_set_num_threads(num_threads);

    double time1 = gbwt::readTimer();

    assert(vg::io::register_libvg_io());

    unique_ptr<handlegraph::HandleGraph> graph = vg::io::VPKG::load_one<handlegraph::HandleGraph>(option_results["graph"].as<string>());
    unique_ptr<gbwt::GBWT> gbwt_index = vg::io::VPKG::load_one<gbwt::GBWT>(option_results["paths"].as<string>());

    PathsIndex paths_index(*gbwt_index, *graph);

    double time2 = gbwt::readTimer();
    cerr << "Load graph and GBWT " << time2 - time1 << " seconds, " << gbwt::inGigabytes(gbwt::memoryUsage()) << " GB" << endl;

    ifstream alignments_istream(option_results["alignments"].as<string>());
    assert(alignments_istream.is_open());

    vector<unordered_map<int32_t, unordered_set<int32_t> > > connected_paths_threads(num_threads);
    vector<vector<vector<AlignmentPath> > > align_paths_threads(num_threads);

    if (is_multipath) {

        AlignmentPathFinder<vg::MultipathAlignment> alignment_path_finder(paths_index, fragment_length_dist.maxLength());

        if (is_single_end) {

             vg::io::for_each_parallel<vg::MultipathAlignment>(alignments_istream, [&](vg::MultipathAlignment & alignment) {

                if (alignment.subpath_size() > 0) {

                    auto align_paths = alignment_path_finder.findAlignmentPaths(alignment);
                    addAlignmentPathsThreaded(&connected_paths_threads, &align_paths_threads, align_paths, omp_get_thread_num());
                }
            });           

        } else {

            vg::io::for_each_interleaved_pair_parallel<vg::MultipathAlignment>(alignments_istream, [&](vg::MultipathAlignment & alignment_1, vg::MultipathAlignment & alignment_2) {

                if (alignment_1.subpath_size() > 0 && alignment_2.subpath_size() > 0) {

                    auto paired_align_paths = alignment_path_finder.findPairedAlignmentPaths(alignment_1, alignment_2);
                    addAlignmentPathsThreaded(&connected_paths_threads, &align_paths_threads, paired_align_paths, omp_get_thread_num());
                }
            });
        }

    } else {

        AlignmentPathFinder<vg::Alignment> alignment_path_finder(paths_index, fragment_length_dist.maxLength());

        if (is_single_end) {

             vg::io::for_each_parallel<vg::Alignment>(alignments_istream, [&](vg::Alignment & alignment) {

                if (alignment.has_path()) {

                    auto align_paths = alignment_path_finder.findAlignmentPaths(alignment);
                    addAlignmentPathsThreaded(&connected_paths_threads, &align_paths_threads, align_paths, omp_get_thread_num());
                }
            });           

        } else {

            vg::io::for_each_interleaved_pair_parallel<vg::Alignment>(alignments_istream, [&](vg::Alignment & alignment_1, vg::Alignment & alignment_2) {

                if (alignment_1.has_path() && alignment_2.has_path()) {

                    auto paired_align_paths = alignment_path_finder.findPairedAlignmentPaths(alignment_1, alignment_2);
                    addAlignmentPathsThreaded(&connected_paths_threads, &align_paths_threads, paired_align_paths, omp_get_thread_num());
                }
            });
        }
    }

    alignments_istream.close();

    double time5 = gbwt::readTimer();
    cerr << "Found alignment paths " << time5 - time2 << " seconds, " << gbwt::inGigabytes(gbwt::memoryUsage()) << " GB" << endl;

    for (size_t i = 1; i < connected_paths_threads.size(); ++i) {

        for (auto & connected_path_clusters: connected_paths_threads.at(i)) {

            auto connected_paths_threads_it = connected_paths_threads.front().emplace(connected_path_clusters.first, unordered_set<int32_t>());
            for (auto & paths: connected_path_clusters.second) {

                connected_paths_threads_it.first->second.emplace(paths);
            }
        }
    }

    double time6 = gbwt::readTimer();
    cerr << "Merged connected path threads " << time6 - time5 << " seconds, " << gbwt::inGigabytes(gbwt::memoryUsage()) << " GB" << endl;

    auto path_clusters = PathClusters(connected_paths_threads.front(), paths_index.index().metadata.haplotype_count);

    double time62 = gbwt::readTimer();
    cerr << "Found path clusters " << time62 - time6 << " seconds, " << gbwt::inGigabytes(gbwt::memoryUsage()) << " GB" << endl;
 
    vector<vector<vector<AlignmentPath> > > clustered_align_paths(path_clusters.cluster_to_path_index.size());

    for (auto & align_paths: align_paths_threads) {

        for (auto & align_path: align_paths) {

            clustered_align_paths.at(path_clusters.path_to_cluster_index.at(align_path.front().ids.front())).push_back(move(align_path));
        }
    }

    double time7 = gbwt::readTimer();
    cerr << "Clustered alignment paths " << time7 - time6 << " seconds, " << gbwt::inGigabytes(gbwt::memoryUsage()) << " GB" << endl;

    vector<vector<ReadPathProbs> > clustered_align_path_probs(path_clusters.cluster_to_path_index.size());
    const double score_log_base = gssw_dna_recover_log_base(1, 4, 0.5, double_precision);

    #pragma omp parallel
    {  
        #pragma omp for
        for (size_t i = 0; i < clustered_align_paths.size(); ++i) {

            clustered_align_path_probs.at(i).reserve(clustered_align_paths.at(i).size());

            unordered_map<int32_t, int32_t> clustered_path_index;
            vector<double> effective_path_lengths; 

            if (!is_long_reads) {

                effective_path_lengths.reserve(path_clusters.cluster_to_path_index.at(i).size());

                for (auto & path_id: path_clusters.cluster_to_path_index.at(i)) {

                    assert(clustered_path_index.emplace(path_id, clustered_path_index.size()).second);
                    effective_path_lengths.emplace_back(paths_index.effectivePathLength(path_id, fragment_length_dist)); 
                }
            }

            for (auto & align_paths: clustered_align_paths.at(i)) {

                clustered_align_path_probs.at(i).emplace_back(ReadPathProbs(clustered_path_index.size(), score_log_base));
                clustered_align_path_probs.at(i).back().calcReadPathProbs(align_paths, clustered_path_index, fragment_length_dist, is_single_end);

                if (!is_long_reads) {

                    clustered_align_path_probs.at(i).back().addPositionalProbs(effective_path_lengths);
                }
            }

            sort(clustered_align_path_probs.at(i).begin(), clustered_align_path_probs.at(i).end());
        }
    }

    double time8 = gbwt::readTimer();
    cerr << "Calculated and sorted alignment path probabilites " << time8 - time7 << " seconds, " << gbwt::inGigabytes(gbwt::memoryUsage()) << " GB" << endl;

    streambuf * output_buffer;
    ofstream output_file;

    if (option_results["output"].as<string>() == "stdout") {

        output_buffer = cout.rdbuf();
    
    } else {

        output_file.open(option_results["output"].as<string>() + ".txt");
        assert(output_file.is_open());

        output_buffer = output_file.rdbuf();
    }

    ostream output_stream(output_buffer);

    for (size_t i = 0; i < clustered_align_path_probs.size(); ++i) {

        const vector<ReadPathProbs> & clustered_probs = clustered_align_path_probs.at(i);

        if (clustered_probs.empty()) {

            continue;
        }

        output_stream << "#\nx Noise";
        for (auto & path_id: path_clusters.cluster_to_path_index.at(i)) {

            output_stream << " " << paths_index.pathName(path_id);
        }
        output_stream << endl;

        output_stream << setprecision(8);

        int32_t num_align_paths = 1;
        int32_t prev_unique_probs_idx = 0;

        for (size_t j = 1; j < clustered_probs.size(); ++j) {

            if (clustered_probs.at(prev_unique_probs_idx) == clustered_probs.at(j)) {

                num_align_paths++;
            
            } else {

                output_stream << num_align_paths << " " << clustered_probs.at(prev_unique_probs_idx) << endl;
                num_align_paths = 1;
                prev_unique_probs_idx = j;
            }
        }

        output_stream << num_align_paths << " " << clustered_probs.at(prev_unique_probs_idx) << endl;
    }

    if (option_results["output"].as<string>() != "stdout") {

        output_file.close();
    }

    double time9 = gbwt::readTimer();
    cerr << "Collapsed and wrote probabilites " << time9 - time8 << " seconds, " << gbwt::inGigabytes(gbwt::memoryUsage()) << " GB" << endl;

	return 0;
}

