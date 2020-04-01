
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
#include "sparsepp/spp.h"
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
#include "read_path_probabilities.hpp"
#include "probability_matrix_writer.hpp"
#include "abundance_estimator.hpp"


static const uint32_t cluster_align_paths_probs_buffer_size = 10;
static const double prob_out_precision = pow(10, -6);

void addAlignmentPathsThreaded(vector<spp::sparse_hash_map<vector<AlignmentPath>, uint32_t> > * all_align_paths_threads, vector<AlignmentPath> * align_paths, const double mean_fragment_length, const uint32_t thread_idx) {

    if (!align_paths->empty()) {

        if (align_paths->size() == 1) {

            align_paths->front().seq_length = mean_fragment_length;
            align_paths->front().score_sum = 1;
        } 

        sort(align_paths->begin(), align_paths->end());

        auto all_align_paths_threads_it = all_align_paths_threads->at(thread_idx).emplace(*align_paths, 0);
        all_align_paths_threads_it.first->second++;
    }    
}

int main(int argc, char* argv[]) {

    cxxopts::Options options("rpvg", "calculate path probabilities and abundances from variation graph read aligments");

    options.add_options("Required")
      ("g,graph", "xg graph file name", cxxopts::value<string>())
      ("p,paths", "GBWT index file name", cxxopts::value<string>())
      ("a,alignments", "gam(p) alignment file name", cxxopts::value<string>())
      ;

    options.add_options("General")
      ("o,output", "output filename", cxxopts::value<string>()->default_value("stdout"))    
      ("t,threads", "number of compute threads", cxxopts::value<uint32_t>()->default_value("1"))
      ("r,seed", "seed for random number generator (default: unix time)", cxxopts::value<int64_t>())
      ("h,help", "print help", cxxopts::value<bool>())
      ;

    options.add_options("Alignment")
      ("u,multipath", "alignment input is multipath gamp format (default: gam)", cxxopts::value<bool>())
      ("s,single-end", "alignment input is single-end reads", cxxopts::value<bool>())
      ("l,long-reads", "alignment input is non-fragmented long reads (single-end only)", cxxopts::value<bool>())
      ;

    options.add_options("Probability")
      ("m,frag-mean", "mean for fragment length distribution", cxxopts::value<double>())
      ("d,frag-sd", "standard deviation for fragment length distribution", cxxopts::value<double>())
      ("b,prob-output", "write read path probabilities to file", cxxopts::value<string>())
      ;

    if (argc == 1) {

        cerr << options.help({"Required", "General", "Alignment", "Probability"}) << endl;
        return 1;
    }

    auto option_results = options.parse(argc, argv);

    if (option_results.count("help")) {

        cerr << options.help({"Required", "General", "Alignment", "Probability"}) << endl;
        return 1;
    }

    if (!option_results.count("graph")) {

        cerr << "ERROR: Graph (xg format) input required (--graph)." << endl;
        return 1;
    }

    if (!option_results.count("paths")) {

        cerr << "ERROR: Paths (GBWT index) input required (--paths)." << endl;
        return 1;
    }

    if (!option_results.count("alignments")) {

        cerr << "ERROR: Alignments (gam or gamp format) input required (--alignments)." << endl;
        return 1;
    }

    uint64_t rng_seed = 0; 

    if (option_results.count("seed")) {

        rng_seed = option_results["seed"].as<uint64_t>();

    } else {

        rng_seed = time(nullptr);
    }

    cerr << "Random number generator seed set to: " << rng_seed << "" << endl;

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

    cerr << endl;

    assert(fragment_length_dist.isValid());

    const uint32_t num_threads = option_results["threads"].as<uint32_t>();

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

    vector<spp::sparse_hash_map<vector<AlignmentPath>, uint32_t> > all_align_paths_threads(num_threads);

    if (is_multipath) {

        AlignmentPathFinder<vg::MultipathAlignment> alignment_path_finder(paths_index, fragment_length_dist.maxLength());

        if (is_single_end) {

             vg::io::for_each_parallel<vg::MultipathAlignment>(alignments_istream, [&](vg::MultipathAlignment & alignment) {

                if (alignment.subpath_size() > 0) {

                    auto align_paths = alignment_path_finder.findAlignmentPaths(alignment);
                    addAlignmentPathsThreaded(&all_align_paths_threads, &align_paths, fragment_length_dist.mean(), omp_get_thread_num());
                }
            });           

        } else {

            vg::io::for_each_interleaved_pair_parallel<vg::MultipathAlignment>(alignments_istream, [&](vg::MultipathAlignment & alignment_1, vg::MultipathAlignment & alignment_2) {

                if (alignment_1.subpath_size() > 0 && alignment_2.subpath_size() > 0) {

                    auto paired_align_paths = alignment_path_finder.findPairedAlignmentPaths(alignment_1, alignment_2);
                    addAlignmentPathsThreaded(&all_align_paths_threads, &paired_align_paths, fragment_length_dist.mean(), omp_get_thread_num());
                }
            });
        }

    } else {

        AlignmentPathFinder<vg::Alignment> alignment_path_finder(paths_index, fragment_length_dist.maxLength());

        if (is_single_end) {

             vg::io::for_each_parallel<vg::Alignment>(alignments_istream, [&](vg::Alignment & alignment) {

                if (alignment.has_path()) {

                    auto align_paths = alignment_path_finder.findAlignmentPaths(alignment);
                    addAlignmentPathsThreaded(&all_align_paths_threads, &align_paths, fragment_length_dist.mean(), omp_get_thread_num());
                }
            });           

        } else {

            vg::io::for_each_interleaved_pair_parallel<vg::Alignment>(alignments_istream, [&](vg::Alignment & alignment_1, vg::Alignment & alignment_2) {

                if (alignment_1.has_path() && alignment_2.has_path()) {

                    auto paired_align_paths = alignment_path_finder.findPairedAlignmentPaths(alignment_1, alignment_2);
                    addAlignmentPathsThreaded(&all_align_paths_threads, &paired_align_paths, fragment_length_dist.mean(), omp_get_thread_num());
                }
            });
        }
    }

    alignments_istream.close();

    double time5 = gbwt::readTimer();
    cerr << "Found alignment paths " << time5 - time2 << " seconds, " << gbwt::inGigabytes(gbwt::memoryUsage()) << " GB" << endl;

    spp::sparse_hash_map<uint32_t, spp::sparse_hash_set<uint32_t> > connected_align_paths(num_threads);

    for (auto & all_align_paths: all_align_paths_threads) {

        cerr << all_align_paths.size() << endl;

        for (auto & align_paths: all_align_paths) {

            auto anchor_path_id = align_paths.first.front().ids.front();

            for (auto & align_path: align_paths.first) {

                for (auto & path_id: align_path.ids) {

                    if (anchor_path_id != path_id) {

                        connected_align_paths[anchor_path_id].emplace(path_id);
                        connected_align_paths[path_id].emplace(anchor_path_id);
                    }
                }
            }
        }
    }

    auto path_clusters = PathClusters(connected_align_paths, paths_index.index().metadata.haplotype_count);

    double time6 = gbwt::readTimer();
    cerr << "Created alignment path cluster index " << time6 - time5 << " seconds, " << gbwt::inGigabytes(gbwt::memoryUsage()) << " GB" << endl;
 
    vector<vector<spp::sparse_hash_map<vector<AlignmentPath>, uint32_t>::iterator> > align_paths_clusters(path_clusters.cluster_to_path_index.size());

    for (size_t i = 0; i < all_align_paths_threads.size(); ++i) {

        auto align_paths_it = all_align_paths_threads.at(i).begin();

        while (align_paths_it != all_align_paths_threads.at(i).end()) {

            align_paths_clusters.at(path_clusters.path_to_cluster_index.at(align_paths_it->first.front().ids.front())).emplace_back(align_paths_it);
            ++align_paths_it;
        }
    }

    double time7 = gbwt::readTimer();
    cerr << "Clustered alignment paths " << time7 - time6 << " seconds, " << gbwt::inGigabytes(gbwt::memoryUsage()) << " GB" << endl;

    vector<vector<vector<pair<ReadPathProbabilities, uint32_t> > > > cluster_align_paths_probs_buffer(num_threads);
    vector<vector<vector<string> > > cluster_path_names_buffer(num_threads);
    vector<vector<vector<double> > > cluster_path_lengths_buffer(num_threads);

    const double score_log_base = gssw_dna_recover_log_base(1, 4, 0.5, double_precision);

    ProbabilityMatrixWriter * prob_matrix_writer = nullptr;
   
    if (option_results.count("prob-output")) {

        prob_matrix_writer = new ProbabilityMatrixWriter(false, option_results["prob-output"].as<string>(), prob_out_precision);
    }

    AbundanceEstimator * em_abundance_estimator = new EMAbundanceEstimator(pow(10, -8), 1000, pow(10, -2));

    #pragma omp parallel
    {  
        #pragma omp for
        for (size_t i = 0; i < align_paths_clusters.size(); ++i) {

            vector<vector<pair<ReadPathProbabilities, uint32_t> > > * thread_cluster_align_paths_probs_buffer = &(cluster_align_paths_probs_buffer.at(omp_get_thread_num()));

            thread_cluster_align_paths_probs_buffer->emplace_back(vector<pair<ReadPathProbabilities, uint32_t> >());            
            thread_cluster_align_paths_probs_buffer->back().reserve(align_paths_clusters.at(i).size());

            unordered_map<uint32_t, uint32_t> clustered_path_index;
            
            vector<vector<double> > * thread_cluster_path_lengths_buffer = &(cluster_path_lengths_buffer.at(omp_get_thread_num()));

            thread_cluster_path_lengths_buffer->emplace_back(vector<double>());            
            thread_cluster_path_lengths_buffer->back().reserve(align_paths_clusters.at(i).size());

            if (is_long_reads) {

                for (auto & path_id: path_clusters.cluster_to_path_index.at(i)) {

                    assert(clustered_path_index.emplace(path_id, clustered_path_index.size()).second);
                    thread_cluster_path_lengths_buffer->back().emplace_back(paths_index.pathLength(path_id)); 
                }
            
            } else {

                for (auto & path_id: path_clusters.cluster_to_path_index.at(i)) {

                    assert(clustered_path_index.emplace(path_id, clustered_path_index.size()).second);
                    thread_cluster_path_lengths_buffer->back().emplace_back(paths_index.effectivePathLength(path_id, fragment_length_dist)); 
                }
            }

            for (auto & align_paths: align_paths_clusters.at(i)) {

                thread_cluster_align_paths_probs_buffer->back().emplace_back(ReadPathProbabilities(clustered_path_index.size(), score_log_base), align_paths->second);
                thread_cluster_align_paths_probs_buffer->back().back().first.calcReadPathProbabilities(align_paths->first, clustered_path_index, fragment_length_dist, is_single_end);

                if (!is_long_reads) {

                    thread_cluster_align_paths_probs_buffer->back().back().first.addPositionalProbabilities(thread_cluster_path_lengths_buffer->back());
                }
            }

            sort(thread_cluster_align_paths_probs_buffer->back().begin(), thread_cluster_align_paths_probs_buffer->back().end());

            vector<vector<string> > * thread_cluster_path_names_buffer = &(cluster_path_names_buffer.at(omp_get_thread_num()));

            thread_cluster_path_names_buffer->emplace_back(vector<string>());            
            thread_cluster_path_names_buffer->back().reserve(align_paths_clusters.at(i).size());

            for (auto & path_id: path_clusters.cluster_to_path_index.at(i)) {

                thread_cluster_path_names_buffer->back().emplace_back(paths_index.pathName(path_id));
            }

            if (thread_cluster_align_paths_probs_buffer->size() == cluster_align_paths_probs_buffer_size) {

                assert(thread_cluster_align_paths_probs_buffer->size() == thread_cluster_path_names_buffer->size());
                assert(thread_cluster_align_paths_probs_buffer->size() == thread_cluster_path_lengths_buffer->size());

                if (prob_matrix_writer) {

                    prob_matrix_writer->lockWriter();

                    for (size_t i = 0; i < cluster_align_paths_probs_buffer_size; ++i) {

                        prob_matrix_writer->writeReadPathProbabilityCluster(thread_cluster_align_paths_probs_buffer->at(i), thread_cluster_path_names_buffer->at(i), thread_cluster_path_lengths_buffer->at(i));
                    }

                    prob_matrix_writer->unlockWriter();
                }

                // for (size_t i = 0; i < cluster_align_paths_probs_buffer_size; ++i) {

                //     em_abundance_estimator->inferClusterAbundance(thread_cluster_align_paths_probs_buffer->at(i), thread_cluster_path_names_buffer->at(i).size());
                // }

                thread_cluster_align_paths_probs_buffer->clear();
                thread_cluster_path_names_buffer->clear();
                thread_cluster_path_lengths_buffer->clear();
            }
        }
    }

    delete prob_matrix_writer;
    delete em_abundance_estimator;

    double time8 = gbwt::readTimer();
    cerr << "Inferred path probabilities and abundances " << time8 - time7 << " seconds, " << gbwt::inGigabytes(gbwt::memoryUsage()) << " GB" << endl;

	return 0;
}

