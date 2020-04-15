
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
#include "path_abundances.hpp"
#include "path_abundance_estimator.hpp"
#include "path_abundance_writer.hpp"


static const uint32_t read_path_cluster_probs_buffer_size = 100;
static const double prob_out_precision = pow(10, -8);

void addAlignmentPathsToIndex(spp::sparse_hash_map<vector<AlignmentPath>, uint32_t> * align_paths_index, vector<AlignmentPath> * align_paths, const double mean_fragment_length) {

    if (!align_paths->empty()) {

        if (align_paths->size() == 1) {

            align_paths->front().seq_length = mean_fragment_length;
            align_paths->front().score_sum = 1;
        } 

        sort(align_paths->begin(), align_paths->end());

        auto threaded_align_paths_index_it = align_paths_index->emplace(*align_paths, 0);
        threaded_align_paths_index_it.first->second++;
    }    
}

spp::sparse_hash_map<string, string> parsePathTranscriptOrigin(const string & filename) {

    spp::sparse_hash_map<string, string> path_transcript_origin;

    ifstream origin_file(filename);
    
    string line;
    string element;

    while (origin_file.good()) {

        getline(origin_file, line);

        if (line.empty()) {

            continue;
        }

        auto line_ss = stringstream(line);

        getline(line_ss, element, '\t');

        if (element == "Name") {

            continue;
        }

        auto path_transcript_origin_it = path_transcript_origin.emplace(element, "");
        assert(path_transcript_origin_it.second);

        getline(line_ss, element, '\t');        
        getline(line_ss, element, '\t');

        path_transcript_origin_it.first->second = element;

        getline(line_ss, element, '\n');
    }

    origin_file.close();

    return path_transcript_origin;
}

int main(int argc, char* argv[]) {

    cxxopts::Options options("fersken", "fersken - infers path probabilities and abundances from variation graph read aligments");

    options.add_options("Required")
      ("g,graph", "xg graph filename", cxxopts::value<string>())
      ("p,paths", "GBWT index filename", cxxopts::value<string>())
      ("a,alignments", "gam(p) alignment filename", cxxopts::value<string>())
      ("i,inference-model", "inference model to use (transcripts, strains or haplotype-transcripts)", cxxopts::value<string>())
      ;

    options.add_options("General")
      ("o,output", "output filename", cxxopts::value<string>()->default_value("stdout"))    
      ("t,threads", "number of compute threads", cxxopts::value<uint32_t>()->default_value("1"))
      ("r,rng-seed", "seed for random number generator (default: unix time)", cxxopts::value<int64_t>())
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

    options.add_options("Abundance")
      ("e,max-em-its", "maximum number of EM iterations", cxxopts::value<uint32_t>()->default_value("1000"))
      ("n,min-abundance", "minimum abundance value", cxxopts::value<double>()->default_value("1e-6"))
      ("y,ploidy", "sample ploidy (used for haplotype-transcript inference, max: 2)", cxxopts::value<uint32_t>()->default_value("2"))
      ("f,num-hap-it", "number of haplotype iterations (used for haplotype-transcript inference)", cxxopts::value<uint32_t>()->default_value("100"))
      ("c,path-origin", "path transcript origin filename (required for haplotype-transcript inference)", cxxopts::value<string>())
      ;

    if (argc == 1) {

        cerr << options.help({"Required", "General", "Alignment", "Probability", "Abundance"}) << endl;
        return 1;
    }

    auto option_results = options.parse(argc, argv);

    if (option_results.count("help")) {

        cerr << options.help({"Required", "General", "Alignment", "Probability", "Abundance"}) << endl;
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

    if (!option_results.count("inference-model")) {

        cerr << "ERROR: Inference model required (--inference-model). Options: transcripts, strains or haplotype-transcripts." << endl;
        return 1;
    }

    const string inference_model = option_results["inference-model"].as<string>();

    if (inference_model != "transcripts" && inference_model != "strains" && inference_model != "haplotype-transcripts") {

        cerr << "ERROR: Inference model provided (--inference-model) not supported. Options: transcripts, strains or haplotype-transcripts." << endl;
        return 1;
    }

    if (option_results["ploidy"].as<uint32_t>() > 2) {

        cerr << "ERROR: Maximum support ploidy (--ploidy) is currently 2." << endl;
        return 1;        
    }

    if (inference_model == "haplotype-transcripts" && !option_results.count("path-origin")) {

        cerr << "ERROR: Path transcript origin information file (--path-origin) needed when running in haplotype-transcript inference mode (--write-info output from vg rna)." << endl;
        return 1;
    }

    uint64_t rng_seed = 0; 

    if (option_results.count("rng-seed")) {

        rng_seed = option_results["rng-seed"].as<uint64_t>();

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

    graph.reset(nullptr);

    double time2 = gbwt::readTimer();
    cerr << "Load graph and GBWT " << time2 - time1 << " seconds, " << gbwt::inGigabytes(gbwt::memoryUsage()) << " GB" << endl;

    ifstream alignments_istream(option_results["alignments"].as<string>());
    assert(alignments_istream.is_open());

    vector<spp::sparse_hash_map<vector<AlignmentPath>, uint32_t> > threaded_align_paths_index(num_threads);

    if (is_multipath) {

        AlignmentPathFinder<vg::MultipathAlignment> alignment_path_finder(paths_index, fragment_length_dist.maxLength());

        if (is_single_end) {

             vg::io::for_each_parallel<vg::MultipathAlignment>(alignments_istream, [&](vg::MultipathAlignment & alignment) {

                if (alignment.subpath_size() > 0) {

                    auto align_paths = alignment_path_finder.findAlignmentPaths(alignment);
                    addAlignmentPathsToIndex(&(threaded_align_paths_index.at(omp_get_thread_num())), &align_paths, fragment_length_dist.mean());
                }
            });           

        } else {

            vg::io::for_each_interleaved_pair_parallel<vg::MultipathAlignment>(alignments_istream, [&](vg::MultipathAlignment & alignment_1, vg::MultipathAlignment & alignment_2) {

                if (alignment_1.subpath_size() > 0 && alignment_2.subpath_size() > 0) {

                    auto paired_align_paths = alignment_path_finder.findPairedAlignmentPaths(alignment_1, alignment_2);
                    addAlignmentPathsToIndex(&(threaded_align_paths_index.at(omp_get_thread_num())), &paired_align_paths, fragment_length_dist.mean());
                }
            });
        }

    } else {

        AlignmentPathFinder<vg::Alignment> alignment_path_finder(paths_index, fragment_length_dist.maxLength());

        if (is_single_end) {

             vg::io::for_each_parallel<vg::Alignment>(alignments_istream, [&](vg::Alignment & alignment) {

                if (alignment.has_path()) {

                    auto align_paths = alignment_path_finder.findAlignmentPaths(alignment);
                    addAlignmentPathsToIndex(&(threaded_align_paths_index.at(omp_get_thread_num())), &align_paths, fragment_length_dist.mean());
                }
            });           

        } else {

            vg::io::for_each_interleaved_pair_parallel<vg::Alignment>(alignments_istream, [&](vg::Alignment & alignment_1, vg::Alignment & alignment_2) {

                if (alignment_1.has_path() && alignment_2.has_path()) {

                    auto paired_align_paths = alignment_path_finder.findPairedAlignmentPaths(alignment_1, alignment_2);
                    addAlignmentPathsToIndex(&(threaded_align_paths_index.at(omp_get_thread_num())), &paired_align_paths, fragment_length_dist.mean());
                }
            });
        }
    }

    alignments_istream.close();

    double time5 = gbwt::readTimer();
    cerr << "Found alignment paths " << time5 - time2 << " seconds, " << gbwt::inGigabytes(gbwt::memoryUsage()) << " GB" << endl;

    spp::sparse_hash_map<uint32_t, spp::sparse_hash_set<uint32_t> > connected_align_paths(num_threads);

    for (auto & align_paths_index: threaded_align_paths_index) {

        cerr << align_paths_index.size() << endl;

        for (auto & align_paths: align_paths_index) {

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

    auto path_clusters = PathClusters(connected_align_paths, paths_index.index().metadata.paths());

    double time6 = gbwt::readTimer();
    cerr << "Created alignment path cluster index " << time6 - time5 << " seconds, " << gbwt::inGigabytes(gbwt::memoryUsage()) << " GB" << endl;
 
    vector<vector<spp::sparse_hash_map<vector<AlignmentPath>, uint32_t>::iterator> > align_paths_clusters(path_clusters.cluster_to_path_index.size());

    for (size_t i = 0; i < threaded_align_paths_index.size(); ++i) {

        auto align_paths_index_it = threaded_align_paths_index.at(i).begin();

        while (align_paths_index_it != threaded_align_paths_index.at(i).end()) {

            align_paths_clusters.at(path_clusters.path_to_cluster_index.at(align_paths_index_it->first.front().ids.front())).emplace_back(align_paths_index_it);
            ++align_paths_index_it;
        }
    }

    double time7 = gbwt::readTimer();
    cerr << "Clustered alignment paths " << time7 - time6 << " seconds, " << gbwt::inGigabytes(gbwt::memoryUsage()) << " GB" << endl;

    vector<vector<vector<pair<ReadPathProbabilities, uint32_t> > > > threaded_read_path_cluster_probs_buffer(num_threads);
    vector<vector<PathAbundances> > threaded_path_cluster_abundances(num_threads);

    for (size_t i = 0; i < num_threads; ++i) {

        threaded_path_cluster_abundances.at(i).reserve(ceil(align_paths_clusters.size()) / static_cast<float>(num_threads));
    }

    ProbabilityMatrixWriter * prob_matrix_writer = nullptr;

    spp::sparse_hash_map<string, string> path_transcript_origin;
   
    if (option_results.count("prob-output")) {

        prob_matrix_writer = new ProbabilityMatrixWriter(false, option_results["prob-output"].as<string>(), prob_out_precision);
    }

    PathAbundanceEstimator * path_abundance_estimator;

    if (inference_model == "transcripts") {

        path_abundance_estimator = new PathAbundanceEstimator(option_results["max-em-its"].as<uint32_t>(), option_results["min-abundance"].as<double>());

    } else if (inference_model == "strains") {

        path_abundance_estimator = new MinimumPathAbundanceEstimator(option_results["max-em-its"].as<uint32_t>(), option_results["min-abundance"].as<double>());

    } else if (inference_model == "haplotype-transcripts") {

        path_abundance_estimator = new NestedPathAbundanceEstimator(option_results["num-hap-it"].as<uint32_t>(), option_results["ploidy"].as<uint32_t>(), rng_seed, option_results["max-em-its"].as<uint32_t>(), option_results["min-abundance"].as<double>());
     
        path_transcript_origin = parsePathTranscriptOrigin(option_results["path-origin"].as<string>());

    } else {

        assert(false);
    }

    const double score_log_base = gssw_dna_recover_log_base(1, 4, 0.5, double_precision);

    auto align_paths_clusters_indices = vector<uint32_t>(align_paths_clusters.size());
    iota(align_paths_clusters_indices.begin(), align_paths_clusters_indices.end(), 0);

    mt19937 mt_rng(rng_seed);
    shuffle(align_paths_clusters_indices.begin(), align_paths_clusters_indices.end(), mt_rng);    

    #pragma omp parallel
    {  
        #pragma omp for
        for (size_t i = 0; i < align_paths_clusters_indices.size(); ++i) {

            auto align_paths_cluster_idx = align_paths_clusters_indices.at(i);

            vector<vector<pair<ReadPathProbabilities, uint32_t> > > * read_path_cluster_probs_buffer = &(threaded_read_path_cluster_probs_buffer.at(omp_get_thread_num()));

            read_path_cluster_probs_buffer->emplace_back(vector<pair<ReadPathProbabilities, uint32_t> >());            
            read_path_cluster_probs_buffer->back().reserve(align_paths_clusters.at(align_paths_cluster_idx).size());

            unordered_map<uint32_t, uint32_t> clustered_path_index;

            vector<Path> cluster_paths;
            cluster_paths.reserve(path_clusters.cluster_to_path_index.at(align_paths_cluster_idx).size());
            
            for (auto & path_id: path_clusters.cluster_to_path_index.at(align_paths_cluster_idx)) {

                assert(clustered_path_index.emplace(path_id, clustered_path_index.size()).second);
                cluster_paths.emplace_back(Path());

                cluster_paths.back().name = paths_index.pathName(path_id);

                auto path_transcript_origin_it = path_transcript_origin.find(cluster_paths.back().name);

                if (path_transcript_origin_it != path_transcript_origin.end()) {

                    cluster_paths.back().origin = path_transcript_origin_it->second;
                }

                cluster_paths.back().length = paths_index.pathLength(path_id); 

                if (is_long_reads) {

                    cluster_paths.back().effective_length = paths_index.pathLength(path_id); 

                } else {

                    cluster_paths.back().effective_length = paths_index.effectivePathLength(path_id, fragment_length_dist); 
                }
            }

            for (auto & align_paths: align_paths_clusters.at(align_paths_cluster_idx)) {

                read_path_cluster_probs_buffer->back().emplace_back(ReadPathProbabilities(clustered_path_index.size(), score_log_base, fragment_length_dist), align_paths->second);
                read_path_cluster_probs_buffer->back().back().first.calcReadPathProbabilities(align_paths->first, clustered_path_index, cluster_paths, is_single_end);
            }

            vector<PathAbundances> * path_cluster_abundances = &(threaded_path_cluster_abundances.at(omp_get_thread_num()));

            path_cluster_abundances->emplace_back(path_abundance_estimator->inferPathClusterAbundances(read_path_cluster_probs_buffer->back(), cluster_paths));

            if (prob_matrix_writer) {

                if (read_path_cluster_probs_buffer->size() == read_path_cluster_probs_buffer_size) {

                    assert(path_cluster_abundances->size() % read_path_cluster_probs_buffer_size == 0);

                    assert(path_cluster_abundances->size() >= read_path_cluster_probs_buffer->size());
                    size_t path_cluster_abundances_idx = path_cluster_abundances->size() - read_path_cluster_probs_buffer->size();

                    prob_matrix_writer->lockWriter();

                    for (size_t j = 0; j < read_path_cluster_probs_buffer->size(); ++j) {

                        prob_matrix_writer->writeReadPathProbabilityCluster(read_path_cluster_probs_buffer->at(j), path_cluster_abundances->at(path_cluster_abundances_idx).paths);
                        ++path_cluster_abundances_idx;
                    }

                    prob_matrix_writer->unlockWriter();
                    read_path_cluster_probs_buffer->clear();       
                } 

            } else {

                read_path_cluster_probs_buffer->clear();
            }
        }
    }

    if (prob_matrix_writer) {

        for (size_t i = 0; i < num_threads; ++i) {

            assert(threaded_path_cluster_abundances.at(i).size() >= threaded_read_path_cluster_probs_buffer.at(i).size());
            size_t path_cluster_abundances_idx = threaded_path_cluster_abundances.at(i).size() - threaded_read_path_cluster_probs_buffer.at(i).size();

            for (size_t j = 0; j < threaded_read_path_cluster_probs_buffer.at(i).size(); ++j) {

                prob_matrix_writer->writeReadPathProbabilityCluster(threaded_read_path_cluster_probs_buffer.at(i).at(j), threaded_path_cluster_abundances.at(i).at(path_cluster_abundances_idx).paths);
                ++path_cluster_abundances_idx;
            }
        }
    } 

    delete prob_matrix_writer;
    delete path_abundance_estimator;

    PathAbundanceWriter path_abundance_writer(option_results["output"].as<string>() == "stdout", option_results["output"].as<string>(), option_results["min-abundance"].as<double>());
    path_abundance_writer.writeThreadedPathClusterAbundances(threaded_path_cluster_abundances);

    double time8 = gbwt::readTimer();
    cerr << "Inferred path probabilities and abundances " << time8 - time7 << " seconds, " << gbwt::inGigabytes(gbwt::memoryUsage()) << " GB" << endl;

	return 0;
}

