
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
#include <thread>
#include <sys/stat.h>

#include "cxxopts.hpp"
#include "gbwt/gbwt.h"
#include "gbwt/fast_locate.h"
#include "sparsepp/spp.h"
#include "vg/io/vpkg.hpp"
#include "vg/io/stream.hpp"
#include "vg/io/basic_stream.hpp"
#include "io/register_libvg_io.hpp"
#include "handlegraph/handle_graph.hpp"
#include "htslib/bgzf.h"
#include "htslib/hts.h"

#include "utils.hpp"
#include "fragment_length_dist.hpp"
#include "paths_index.hpp"
#include "alignment_path.hpp"
#include "alignment_path_finder.hpp"
#include "producer_consumer_queue.hpp"
#include "path_clusters.hpp"
#include "read_path_probabilities.hpp"
#include "path_estimator.hpp"
#include "path_posterior_estimator.hpp"
#include "path_abundance_estimator.hpp"
#include "path_cluster_estimates.hpp"
#include "threaded_output_writer.hpp"

const uint32_t align_paths_buffer_size = 10000;
const uint32_t fragment_length_min_mapq = 40;

typedef spp::sparse_hash_map<vector<AlignmentPath>, uint32_t> align_paths_index_t;
typedef spp::sparse_hash_map<uint32_t, spp::sparse_hash_set<uint32_t> > connected_align_paths_t;

typedef ProducerConsumerQueue<vector<vector<AlignmentPath> > *> align_paths_buffer_queue_t;


void addAlignmentPathsToBuffer(const vector<AlignmentPath> & align_paths, vector<vector<AlignmentPath> > * align_paths_buffer) {

    if (!align_paths.empty()) {

        assert(align_paths.size() > 1);

        if (align_paths.size() == 2) {

            align_paths_buffer->emplace_back(align_paths);
        
        } else {

            auto align_paths_it = align_paths.begin();
            assert(align_paths_it != align_paths.end());

            align_paths_buffer->emplace_back(vector<AlignmentPath>(1, *align_paths_it));
            ++align_paths_it;

            while (align_paths_it != align_paths.end()) {

                const AlignmentPath & last_align_paths_buffer = align_paths_buffer->back().back();

                if (last_align_paths_buffer.gbwt_search == align_paths_it->gbwt_search && last_align_paths_buffer.is_multimap == align_paths_it->is_multimap && last_align_paths_buffer.frag_length == align_paths_it->frag_length && last_align_paths_buffer.min_mapq == align_paths_it->min_mapq) {

                    assert(last_align_paths_buffer.score_sum >= align_paths_it->score_sum);

                } else {

                    align_paths_buffer->back().emplace_back(*align_paths_it);
                }

                ++align_paths_it;
            }

            assert(align_paths_buffer->back().size() > 1);
        }
    }
}

template<class AlignmentType> 
void findAlignmentPaths(ifstream & alignments_istream, align_paths_buffer_queue_t * align_paths_buffer_queue, const AlignmentPathFinder<AlignmentType> & align_path_finder, const uint32_t num_threads) {

    auto threaded_align_paths_buffer = vector<vector<vector<AlignmentPath > > *>(num_threads);

    for (auto & align_paths_buffer: threaded_align_paths_buffer) {

        align_paths_buffer = new vector<vector<AlignmentPath > >();
        align_paths_buffer->reserve(align_paths_buffer_size);
    }
  
    vg::io::for_each_parallel<AlignmentType>(alignments_istream, [&](AlignmentType & alignment) {

        vector<vector<AlignmentPath > > * align_paths_buffer = threaded_align_paths_buffer.at(omp_get_thread_num());
        addAlignmentPathsToBuffer(align_path_finder.findAlignmentPaths(alignment), align_paths_buffer);

        if (align_paths_buffer->size() == align_paths_buffer_size) {

            align_paths_buffer_queue->push(align_paths_buffer);
            
            threaded_align_paths_buffer.at(omp_get_thread_num()) = new vector<vector<AlignmentPath > >();
            threaded_align_paths_buffer.at(omp_get_thread_num())->reserve(align_paths_buffer_size);
        }
    });

    for (auto & align_paths_buffer: threaded_align_paths_buffer) {

        align_paths_buffer_queue->push(align_paths_buffer);
    }
}

template<class AlignmentType> 
void findPairedAlignmentPaths(ifstream & alignments_istream, align_paths_buffer_queue_t * align_paths_buffer_queue, const AlignmentPathFinder<AlignmentType> & align_path_finder, const uint32_t num_threads) {

    auto threaded_align_paths_buffer = vector<vector<vector<AlignmentPath > > *>(num_threads);

    for (auto & align_paths_buffer: threaded_align_paths_buffer) {

        align_paths_buffer = new vector<vector<AlignmentPath > >();
        align_paths_buffer->reserve(align_paths_buffer_size);
    }
  
    vg::io::for_each_interleaved_pair_parallel<AlignmentType>(alignments_istream, [&](AlignmentType & alignment_1, AlignmentType & alignment_2) {

        vector<vector<AlignmentPath > > * align_paths_buffer = threaded_align_paths_buffer.at(omp_get_thread_num());
        addAlignmentPathsToBuffer(align_path_finder.findPairedAlignmentPaths(alignment_1, alignment_2), align_paths_buffer);

        if (align_paths_buffer->size() == align_paths_buffer_size) {

            align_paths_buffer_queue->push(align_paths_buffer);

            threaded_align_paths_buffer.at(omp_get_thread_num()) = new vector<vector<AlignmentPath > >();
            threaded_align_paths_buffer.at(omp_get_thread_num())->reserve(align_paths_buffer_size);
        }
    });

    for (auto & align_paths_buffer: threaded_align_paths_buffer) {

        align_paths_buffer_queue->push(align_paths_buffer);
    }
}

void addAlignmentPathsBufferToIndexes(align_paths_buffer_queue_t * align_paths_buffer_queue, align_paths_index_t * align_paths_index, FragmentLengthDist * fragment_length_dist, const uint32_t mean_pre_fragment_length) {

    vector<vector<AlignmentPath> > * align_paths_buffer = nullptr;
    vector<uint32_t> fragment_length_counts(1000, 0);

    while (align_paths_buffer_queue->pop(&align_paths_buffer)) {

        for (auto & align_paths: *align_paths_buffer) {

            assert(align_paths.size() > 1);
            assert(align_paths.back().frag_length == 0);

            if (align_paths.front().min_mapq >= fragment_length_min_mapq && !align_paths.front().is_multimap) {

                uint32_t cur_fragment_length = align_paths.front().frag_length;
                bool cur_length_is_constant = true;

                for (size_t j = 1; j < align_paths.size() - 1; ++j) {

                    assert(align_paths.at(j).min_mapq >= fragment_length_min_mapq);
                    assert(!align_paths.at(j).is_multimap);

                    if (align_paths.at(j).frag_length != cur_fragment_length) {

                        cur_length_is_constant = false;
                        break;
                    }
                }

                if (cur_length_is_constant) {

                    if (fragment_length_counts.size() <= cur_fragment_length) {
                        
                        fragment_length_counts.resize(cur_fragment_length + 1, 0);
                    }

                    fragment_length_counts.at(cur_fragment_length)++;
                }   
            }

            if (align_paths.size() == 2) {       

                align_paths.front().frag_length = mean_pre_fragment_length;      
                align_paths.front().score_sum = 1;       
            } 

            auto threaded_align_paths_index_it = align_paths_index->emplace(align_paths, 0);
            threaded_align_paths_index_it.first->second++;
        } 

        delete align_paths_buffer;
    }

    *fragment_length_dist = FragmentLengthDist(fragment_length_counts);
}

spp::sparse_hash_map<string, PathInfo> parseHaplotypeTranscriptInfo(const string & filename, const bool parse_haplotype_ids) {

    spp::sparse_hash_map<string, PathInfo> haplotype_transcript_info;

    ifstream info_file(filename);
    
    string line;
    string element;

    spp::sparse_hash_map<string, uint32_t> transcript_id_index;
    spp::sparse_hash_map<string, uint32_t> haplotype_id_index;

    while (info_file.good()) {

        getline(info_file, line);

        if (line.empty()) {

            continue;
        }

        auto line_ss = stringstream(line);

        getline(line_ss, element, '\t');

        if (element == "Name") {

            continue;
        }

        auto haplotype_transcript_info_it = haplotype_transcript_info.emplace(element, PathInfo(element));
        assert(haplotype_transcript_info_it.second);

        getline(line_ss, element, '\t');        
        getline(line_ss, element, '\t');

        auto transcript_id_index_it = transcript_id_index.emplace(element, transcript_id_index.size());
        haplotype_transcript_info_it.first->second.group_id = transcript_id_index_it.first->second;

        getline(line_ss, element, '\t');
        getline(line_ss, element, '\n');

        if (parse_haplotype_ids) {

            for (auto & haplotype: Utils::splitString(element, ',')) {

                auto haplotype_id_index_it = haplotype_id_index.emplace(haplotype, haplotype_id_index.size());
                assert(haplotype_transcript_info_it.first->second.source_ids.emplace(haplotype_id_index_it.first->second).second);
            }

            haplotype_transcript_info_it.first->second.source_count = haplotype_transcript_info_it.first->second.source_ids.size();

        } else {

            haplotype_transcript_info_it.first->second.source_count = count(element.begin(), element.end(), ',') + 1;
        }
    }

    info_file.close();
    return haplotype_transcript_info;
}

// https://stackoverflow.com/questions/12774207/fastest-way-to-check-if-a-file-exist-using-standard-c-c11-c
bool doesFileExist(const string name) {
  
    struct stat buffer;
    return (stat(name.c_str(), &buffer) == 0); 
}

int main(int argc, char* argv[]) {

    cxxopts::Options options("rpvg", "rpvg - infers path posterior probabilities and abundances from variation graph read alignments");

    options.add_options("Required")
      ("g,graph", "xg graph filename", cxxopts::value<string>())
      ("p,paths", "GBWT index filename", cxxopts::value<string>())
      ("a,alignments", "gam(p) alignment filename", cxxopts::value<string>())
      ("o,output-prefix", "prefix used for output filenames (e.g. <prefix>.txt)", cxxopts::value<string>())
      ("i,inference-model", "inference model to use (haplotypes, transcripts, strains or haplotype-transcripts)", cxxopts::value<string>())
      ;

    options.add_options("General")
      ("t,threads", "number of compute threads (+= 1 I/O thread)", cxxopts::value<uint32_t>()->default_value("1"))
      ("r,rng-seed", "seed for random number generator (default: unix time)", cxxopts::value<uint64_t>())
      ("h,help", "print help", cxxopts::value<bool>())
      ;

    options.add_options("Alignment")
      ("e,strand-specific", "strand-specific library type (fr: read1 forward, rf: read1 reverse)", cxxopts::value<string>()->default_value("unstranded"))
      ("u,single-path", "alignment input is single-path gam format (default: multipath gamp)", cxxopts::value<bool>())
      ("s,single-end", "alignment input is single-end reads", cxxopts::value<bool>())
      ("l,long-reads", "alignment input is single-molecule long reads (single-end only)", cxxopts::value<bool>())
      ;

    options.add_options("Probability")
      ("m,frag-mean", "mean for fragment length distribution", cxxopts::value<double>())
      ("d,frag-sd", "standard deviation for fragment length distribution", cxxopts::value<double>())
      ("b,write-probs", "write read path probabilities to file (<prefix>_probs.txt.gz)", cxxopts::value<bool>())
      ("max-par-offset", "maximum start and end offset allowed for partial path alignments", cxxopts::value<uint32_t>()->default_value("4"))
      ("est-missing-prob", "estimate the probability that the correct alignment path is missing (experimental)", cxxopts::value<bool>())
      ("filt-best-score", "filter alignments with a best score fraction of <value> below optimal", cxxopts::value<double>()->default_value("0.9"))
      ("min-noise-prob", "minimum probability that alignment is incorrect", cxxopts::value<double>()->default_value("1e-4"))
      ("prob-precision", "precision threshold used to collapse similar probabilities and filter output", cxxopts::value<double>()->default_value("1e-8"))
      ("path-node-cluster", "also cluster paths sharing a node (default: paths sharing a read)", cxxopts::value<bool>())
      ;

    options.add_options("Haplotyping")
      ("y,ploidy", "max sample ploidy", cxxopts::value<uint32_t>()->default_value("2"))
      ("f,path-info", "path haplotype/transcript info filename (required for haplotype-transcript inference)", cxxopts::value<string>())
      ("ind-hap-inference", "infer haplotypes independently for each transcript in haplotype-transcript inference", cxxopts::value<bool>())
      ("num-hap-samples", "number of haplotyping samples in haplotype-transcript inference", cxxopts::value<uint32_t>()->default_value("1000"))
      ("use-hap-gibbs", "use Gibbs sampling for haplotype inference", cxxopts::value<bool>())
      ;

    options.add_options("Quantification")
      ("n,num-gibbs-samples", "number of Gibbs samples per haplotype sample (written to <prefix>_gibbs.txt.gz)", cxxopts::value<uint32_t>()->default_value("0"))
      ("max-em-its", "maximum number of quantification EM iterations", cxxopts::value<uint32_t>()->default_value("10000"))
      ("max-rel-em-conv", "maximum relative abundance used for EM convergence", cxxopts::value<double>()->default_value("0.01"))
      ("gibbs-thin-its", "number of Gibbs iterations between samples", cxxopts::value<uint32_t>()->default_value("25"))      
      ;

    if (argc == 1) {

        cerr << options.help({"Required", "General", "Alignment", "Probability", "Haplotyping", "Quantification"}) << endl;
        return 1;
    }

    auto option_results = options.parse(argc, argv);

    if (option_results.count("help")) {

        cerr << options.help({"Required", "General", "Alignment", "Probability", "Haplotyping", "Quantification"}) << endl;
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

    if (!option_results.count("output-prefix")) {

        cerr << "ERROR: Prefix used for output filenames required (--output-prefix)." << endl;
        return 1;
    }

    if (!option_results.count("inference-model")) {

        cerr << "ERROR: Inference model required (--inference-model). Options: haplotypes, transcripts, strains or haplotype-transcripts." << endl;
        return 1;
    }

    const string inference_model = option_results["inference-model"].as<string>();

    if (inference_model != "haplotypes" && inference_model != "transcripts" && inference_model != "strains" && inference_model != "haplotype-transcripts") {

        cerr << "ERROR: Inference model provided (--inference-model) not supported. Options: haplotypes, transcripts, strains or haplotype-transcripts." << endl;
        return 1;
    }

    const uint32_t num_threads = option_results["threads"].as<uint32_t>();
    assert(num_threads > 0);

    omp_set_num_threads(num_threads);

    uint64_t rng_seed = 0; 

    if (option_results.count("rng-seed")) {

        rng_seed = option_results["rng-seed"].as<uint64_t>();

    } else {

        rng_seed = time(nullptr);
    }


    cerr << "Running rpvg (commit: " << GIT_COMMIT << ")" << endl;
    cerr << "Random number generator seed: " << rng_seed << endl;

    const string library_type = option_results["strand-specific"].as<string>();

    if (library_type != "unstranded" && library_type != "fr" && library_type != "rf") {

        cerr << "ERROR: Strand-specific library type provided (--strand-specific) not supported. Options: unstranded, fr or rf." << endl;
        return 1;
    }

    const bool is_single_end = (option_results.count("single-end") || option_results.count("long-reads"));
    const bool is_long_reads = option_results.count("long-reads");
    const bool is_single_path = option_results.count("single-path");

    if (option_results.count("frag-mean") != option_results.count("frag-sd")) {

        cerr << "ERROR: Both --frag-mean and --frag-sd needs to be given as input. Alternative, no values can be given for paired-end, non-long read alignments and the parameter estimated during mapping will be used instead (contained in the alignment file)." << endl;
        return 1;
    }

    FragmentLengthDist pre_fragment_length_dist; 

    if (is_long_reads) {

        assert(is_single_end);
        pre_fragment_length_dist = FragmentLengthDist(1, 1);

    } else if (!option_results.count("frag-mean") && !option_results.count("frag-sd")) {

        if (is_single_end) {

            cerr << "ERROR: Both --frag-mean and --frag-sd needs to be given as input when using single-end, short read alignments." << endl;
            return 1;            
        }

        ifstream frag_alignments_istream(option_results["alignments"].as<string>());
        assert(frag_alignments_istream.is_open());

        pre_fragment_length_dist = FragmentLengthDist(&frag_alignments_istream, !is_single_path);

        frag_alignments_istream.close();

        if (!pre_fragment_length_dist.isValid()) {

            cerr << "ERROR: No fragment length distribution parameters found in alignments. Use --frag-mean and --frag-sd instead." << endl;
            return 1;
        
        } else {

            cerr << "Fragment length distribution parameters found in alignment (mean: " << pre_fragment_length_dist.mean() << ", standard deviation: " << pre_fragment_length_dist.sd() << ")" << endl;
        }      

    } else {

        pre_fragment_length_dist = FragmentLengthDist(option_results["frag-mean"].as<double>(), option_results["frag-sd"].as<double>());

        cerr << "Fragment length distribution parameters given as input (mean: " << pre_fragment_length_dist.mean() << ", standard deviation: " << pre_fragment_length_dist.sd() << ")" << endl;
    }

    cerr << endl;

    const uint32_t max_partial_offset = option_results["max-par-offset"].as<uint32_t>();
    const bool est_missing_noise_prob = option_results.count("est-missing-prob");
  
    const double min_best_score_filter = option_results["filt-best-score"].as<double>();
    assert(min_best_score_filter >= 0 && min_best_score_filter <= 1);

    const double min_noise_prob = option_results["min-noise-prob"].as<double>();
    assert(min_noise_prob >= 0 && min_noise_prob <= 1);

    const double prob_precision = option_results["prob-precision"].as<double>();
    assert(prob_precision >= 0 && prob_precision <= 1);

    const uint32_t ploidy = option_results["ploidy"].as<uint32_t>();

    if (ploidy == 0) {

        cerr << "ERROR: Ploidy (--ploidy) can not be 0." << endl;
        return 1;        
    }

    if (inference_model == "haplotype-transcripts" && !option_results.count("path-info")) {

        cerr << "ERROR: Path haplotype/transcript information file (--path-info) needed when running in haplotype-transcript inference mode (--write-info output from vg rna)." << endl;
        return 1;
    }

    assert(pre_fragment_length_dist.isValid());

    const bool ind_hap_inference = option_results.count("ind-hap-inference");
    const uint32_t num_hap_samples = option_results["num-hap-samples"].as<uint32_t>();
    const bool use_hap_gibbs = option_results.count("use-hap-gibbs");

    const uint32_t num_gibbs_samples = option_results["num-gibbs-samples"].as<uint32_t>();
    const uint32_t max_em_its = option_results["max-em-its"].as<uint32_t>();
    const double max_rel_em_conv = option_results["max-rel-em-conv"].as<double>();
    const uint32_t gibbs_thin_its = option_results["gibbs-thin-its"].as<uint32_t>();

    double time_init = gbwt::readTimer();

    assert(vg::io::register_libvg_io());

    unique_ptr<handlegraph::HandleGraph> graph = vg::io::VPKG::load_one<handlegraph::HandleGraph>(option_results["graph"].as<string>());
    unique_ptr<gbwt::GBWT> gbwt_index = vg::io::VPKG::load_one<gbwt::GBWT>(option_results["paths"].as<string>());

    unique_ptr<gbwt::FastLocate> r_index;

    if (doesFileExist(option_results["paths"].as<string>() + ".ri")) {

        r_index = move(vg::io::VPKG::load_one<gbwt::FastLocate>(option_results["paths"].as<string>() + ".ri"));
        r_index->setGBWT(*gbwt_index);        

    } else {

        r_index = std::make_unique<gbwt::FastLocate>();
    }

    PathsIndex paths_index(*gbwt_index, *r_index, *graph);
    graph.reset(nullptr);

    if (paths_index.numberOfPaths() == 0) {

        cerr << "ERROR: The GBWT index does not contain any paths." << endl;
        return 1;        
    }

    double time_load = gbwt::readTimer();

    if (r_index->empty()) {

        cerr << "Loaded graph and GBWT (" << time_load - time_init << " seconds, " << gbwt::inGigabytes(gbwt::memoryUsage()) << " GB)" << endl;

    } else {

        cerr << "Loaded graph, GBWT and r-index (" << time_load - time_init << " seconds, " << gbwt::inGigabytes(gbwt::memoryUsage()) << " GB)" << endl;        
    }

    ifstream alignments_istream(option_results["alignments"].as<string>());
    assert(alignments_istream.is_open());

    align_paths_index_t align_paths_index;
    auto align_paths_buffer_queue = new align_paths_buffer_queue_t(num_threads * 3);

    FragmentLengthDist fragment_length_dist;

    thread indexing_thread(addAlignmentPathsBufferToIndexes, align_paths_buffer_queue, &align_paths_index, &fragment_length_dist, pre_fragment_length_dist.mean());

    if (is_single_path) {
        
        AlignmentPathFinder<vg::Alignment> align_path_finder(paths_index, library_type, pre_fragment_length_dist.maxLength(), max_partial_offset, est_missing_noise_prob, min_best_score_filter);

        if (is_single_end) {

            findAlignmentPaths<vg::Alignment>(alignments_istream, align_paths_buffer_queue, align_path_finder, num_threads);

        } else {

            findPairedAlignmentPaths<vg::Alignment>(alignments_istream, align_paths_buffer_queue, align_path_finder, num_threads);
        }

    } else {

        AlignmentPathFinder<vg::MultipathAlignment> align_path_finder(paths_index, library_type, pre_fragment_length_dist.maxLength(), max_partial_offset, est_missing_noise_prob, min_best_score_filter);

        if (is_single_end) {

            findAlignmentPaths<vg::MultipathAlignment>(alignments_istream, align_paths_buffer_queue, align_path_finder, num_threads);

        } else {

            findPairedAlignmentPaths<vg::MultipathAlignment>(alignments_istream, align_paths_buffer_queue, align_path_finder, num_threads);
        }        
    }

    alignments_istream.close();
    align_paths_buffer_queue->pushedLast();

    indexing_thread.join();
    delete align_paths_buffer_queue;

    cerr << align_paths_index.size() << endl;

    if (is_single_end || is_long_reads) {

        fragment_length_dist = pre_fragment_length_dist;

    } else {

        if (!fragment_length_dist.isValid()) {

            if (option_results.count("frag-mean") && option_results.count("frag-sd")) {

                cerr << "Warning: Less than 2 unambiguous read pairs available to re-estimate fragment length distribution parameters from alignment paths. Will use parameters given as input instead (mean: " << pre_fragment_length_dist.mean() << ", standard deviation: " << pre_fragment_length_dist.sd() << ")" << endl;

                fragment_length_dist = pre_fragment_length_dist;

            } else {

                cerr << "Error: Less than 2 unambiguous read pairs available to re-estimate fragment length distribution parameters from alignment paths. Use --frag-mean and --frag-sd instead." << endl;
                return 1;
            }
        
        } else {

            cerr << "Fragment length distribution parameters re-estimated from alignment paths (mean: " << fragment_length_dist.mean() << ", standard deviation: " << fragment_length_dist.sd() << ")" << endl;
        }
    }

    double time_align = gbwt::readTimer();
    cerr << "Found alignment paths (" << time_align - time_load << " seconds, " << gbwt::inGigabytes(gbwt::memoryUsage()) << " GB)" << endl;

    PathClusters path_clusters(num_threads, paths_index, align_paths_index);

    if (option_results.count("path-node-cluster")) {

        path_clusters.addNodeClusters(paths_index);
    }

    cerr << path_clusters.cluster_to_paths_index.size() << endl;

    vector<vector<vector<align_paths_index_t::iterator> > > align_paths_clusters(path_clusters.cluster_to_paths_index.size(), vector<vector<align_paths_index_t::iterator> >(num_threads));

    #pragma omp parallel num_threads(num_threads)
    {
        #pragma omp for schedule(static, 1)
        for (size_t i = 0; i < num_threads; ++i) {

            auto align_paths_index_it = align_paths_index.begin();
            advance(align_paths_index_it, i);

            uint32_t cur_align_paths_idx = i; 

            while (cur_align_paths_idx < align_paths_index.size()) {

                assert(!align_paths_index_it->first.front().gbwt_search.first.empty());
                const uint32_t anchor_path_id = paths_index.locatePathIds(align_paths_index_it->first.front().gbwt_search).front();
                
                align_paths_clusters.at(path_clusters.path_to_cluster_index.at(anchor_path_id)).at(i).emplace_back(align_paths_index_it);

                advance(align_paths_index_it, num_threads);
                cur_align_paths_idx += num_threads;
            }
        }
    }

    double time_clust = gbwt::readTimer();
    cerr << "Clustered alignment paths (" << time_clust - time_align << " seconds, " << gbwt::inGigabytes(gbwt::memoryUsage()) << " GB)" << endl;

    spp::sparse_hash_map<string, PathInfo> haplotype_transcript_info;

    PathEstimator * path_estimator;

    if (inference_model == "haplotypes") {

        path_estimator = new PathGroupPosteriorEstimator(ploidy, use_hap_gibbs, prob_precision);

    } else if (inference_model == "transcripts") {

        path_estimator = new PathAbundanceEstimator(max_em_its, max_rel_em_conv, num_gibbs_samples, gibbs_thin_its, prob_precision);

    } else if (inference_model == "strains") {

        path_estimator = new MinimumPathAbundanceEstimator(max_em_its, max_rel_em_conv, num_gibbs_samples, gibbs_thin_its, prob_precision);

    } else if (inference_model == "haplotype-transcripts") {

        path_estimator = new NestedPathAbundanceEstimator(ploidy, num_hap_samples, !ind_hap_inference, use_hap_gibbs, max_em_its, max_rel_em_conv, num_gibbs_samples, gibbs_thin_its, prob_precision);
        haplotype_transcript_info = parseHaplotypeTranscriptInfo(option_results["path-info"].as<string>(), !ind_hap_inference);

    } else {

        assert(false);
    }

    ProbabilityClusterWriter * prob_cluster_writer = nullptr;

    if (option_results.count("write-probs")) {

        prob_cluster_writer = new ProbabilityClusterWriter(option_results["output-prefix"].as<string>(), num_threads, prob_precision);
    }

    GibbsSamplesWriter * gibbs_samples_writer = nullptr;

    if (num_gibbs_samples > 0) {

        gibbs_samples_writer = new GibbsSamplesWriter(option_results["output-prefix"].as<string>(), num_threads, num_gibbs_samples);
    }

    vector<vector<pair<uint32_t, PathClusterEstimates> > > threaded_path_cluster_estimates(num_threads);

    for (size_t i = 0; i < num_threads; ++i) {

        threaded_path_cluster_estimates.at(i).reserve(ceil(align_paths_clusters.size()) / static_cast<float>(num_threads));
    }

    auto align_paths_clusters_indices = vector<pair<uint32_t, uint32_t> >();
    align_paths_clusters_indices.reserve(align_paths_clusters.size());

    for (size_t i = 0; i < align_paths_clusters.size(); ++i) {

        uint32_t num_align_paths = 0;
        assert(align_paths_clusters.at(i).size() ==  num_threads);

        for (auto & align_paths: align_paths_clusters.at(i)) {

            num_align_paths += align_paths.size();
        }

        align_paths_clusters_indices.emplace_back(num_align_paths, i);
    }

    sort(align_paths_clusters_indices.rbegin(), align_paths_clusters_indices.rend());

    #pragma omp parallel for schedule(dynamic, 1)
    for (size_t i = 0; i < align_paths_clusters_indices.size(); ++i) {

        auto align_paths_cluster_idx = align_paths_clusters_indices.at(i).second;
        auto thread_id = omp_get_thread_num();

        // double debug_time = gbwt::readTimer();

        // if (path_clusters.cluster_to_paths_index.at(align_paths_cluster_idx).size() > 1000 || align_paths_clusters.at(align_paths_cluster_idx).size() > 1000) {

        //     #pragma omp critical
        //     {
                
        //         cerr << "DEBUG: Start " << omp_get_thread_num() << ": " << i << " " << path_clusters.cluster_to_paths_index.at(align_paths_cluster_idx).size() << " " << align_paths_clusters.at(align_paths_cluster_idx).size() << " " << gbwt::inGigabytes(gbwt::memoryUsage()) << endl;
        //     }
        // }

        spp::sparse_hash_map<uint32_t, uint32_t> clustered_path_index;

        auto * path_cluster_estimates = &(threaded_path_cluster_estimates.at(omp_get_thread_num()));
        path_cluster_estimates->emplace_back(i + 1, PathClusterEstimates());

        path_cluster_estimates->back().second.paths.reserve(path_clusters.cluster_to_paths_index.at(align_paths_cluster_idx).size());
        
        for (auto & path_id: path_clusters.cluster_to_paths_index.at(align_paths_cluster_idx)) {

            assert(clustered_path_index.emplace(path_id, clustered_path_index.size()).second);

            if (inference_model == "haplotype-transcripts") {

                auto haplotype_transcript_info_it = haplotype_transcript_info.find(paths_index.pathName(path_id));
                assert(haplotype_transcript_info_it != haplotype_transcript_info.end());

                path_cluster_estimates->back().second.paths.emplace_back(move(haplotype_transcript_info_it->second));
            
            } else {

                path_cluster_estimates->back().second.paths.emplace_back(PathInfo(paths_index.pathName(path_id)));
            }

            path_cluster_estimates->back().second.paths.back().length = paths_index.pathLength(path_id); 

            if (is_long_reads) {

                path_cluster_estimates->back().second.paths.back().effective_length = paths_index.pathLength(path_id); 

            } else {

                path_cluster_estimates->back().second.paths.back().effective_length = paths_index.effectivePathLength(path_id, fragment_length_dist); 
            }
        }

        vector<ReadPathProbabilities> read_path_cluster_probs;
        read_path_cluster_probs.reserve(align_paths_clusters_indices.at(i).first);

        for (auto & threaded_align_paths: align_paths_clusters.at(align_paths_cluster_idx)) {

            for (auto & align_paths: threaded_align_paths) {

                vector<vector<gbwt::size_type> > align_paths_ids;
                align_paths_ids.reserve(align_paths->first.size());

                for (auto & align_path: align_paths->first) {

                    align_paths_ids.emplace_back(paths_index.locatePathIds(align_path.gbwt_search));
                }

                read_path_cluster_probs.emplace_back(ReadPathProbabilities(align_paths->second, prob_precision));
                read_path_cluster_probs.back().calcAlignPathProbs(align_paths->first, align_paths_ids, clustered_path_index, path_cluster_estimates->back().second.paths, fragment_length_dist, is_single_end, min_noise_prob);
            }
        }

        sort(read_path_cluster_probs.begin(), read_path_cluster_probs.end());

        if (!read_path_cluster_probs.empty()) {        

            uint32_t prev_unique_probs_idx = 0;

            for (size_t i = 1; i < read_path_cluster_probs.size(); ++i) {

                if (!read_path_cluster_probs.at(prev_unique_probs_idx).quickMergeIdentical(read_path_cluster_probs.at(i))) {

                    if (prev_unique_probs_idx + 1 < i) {

                        read_path_cluster_probs.at(prev_unique_probs_idx + 1) = read_path_cluster_probs.at(i);
                    }

                    prev_unique_probs_idx++;
                }
            }

            read_path_cluster_probs.resize(prev_unique_probs_idx + 1);
        }

        // Need better solution for this
        mt19937 mt_rng = mt19937(rng_seed + i);
        path_estimator->estimate(&(path_cluster_estimates->back().second), read_path_cluster_probs, &mt_rng);

        if (prob_cluster_writer) {

            prob_cluster_writer->addCluster(read_path_cluster_probs, path_cluster_estimates->back().second.paths);
        } 

        if (gibbs_samples_writer) {

            gibbs_samples_writer->addSamples(path_cluster_estimates->back());
            path_cluster_estimates->back().second.gibbs_read_count_samples.clear();
        }

        // if (path_clusters.cluster_to_paths_index.at(align_paths_cluster_idx).size() > 1000 || align_paths_clusters.at(align_paths_cluster_idx).size() > 1000) {

        //     #pragma omp critical
        //     {
                
        //         cerr << "DEBUG: End " << omp_get_thread_num() << ": " << i << " " << path_clusters.cluster_to_paths_index.at(align_paths_cluster_idx).size() << " " << align_paths_clusters.at(align_paths_cluster_idx).size() << " " << gbwt::inGigabytes(gbwt::memoryUsage()) << " " << gbwt::readTimer() - debug_time << endl;
        //     }
        // }
    }

    delete path_estimator;

    if (prob_cluster_writer) {

        prob_cluster_writer->close();
    } 

    if (gibbs_samples_writer) {

        gibbs_samples_writer->close();
    }

    delete prob_cluster_writer;
    delete gibbs_samples_writer;

    if (inference_model == "haplotypes") {

        PosteriorEstimatesWriter posterior_estimates_writer(option_results["output-prefix"].as<string>(), num_threads, ploidy, prob_precision);

        for (auto & path_cluster_estimates: threaded_path_cluster_estimates) {

            posterior_estimates_writer.addEstimates(path_cluster_estimates);
        }

        posterior_estimates_writer.close();

    } else {

        double total_transcript_count = 0;

        for (auto & path_cluster_estimates_thread: threaded_path_cluster_estimates) {

            for (auto & path_cluster_estimates: path_cluster_estimates_thread) {

                assert(path_cluster_estimates.second.paths.size() == path_cluster_estimates.second.abundances.cols());

                for (size_t i = 0; i < path_cluster_estimates.second.paths.size(); ++i) {

                    if (path_cluster_estimates.second.paths.at(i).effective_length > 0) {

                        total_transcript_count += (path_cluster_estimates.second.abundances(0, i) / path_cluster_estimates.second.paths.at(i).effective_length);
                    }
                }
            }
        }

        AbundanceEstimatesWriter abundance_estimates_writer(option_results["output-prefix"].as<string>(), num_threads, total_transcript_count);

        for (auto & path_cluster_estimates: threaded_path_cluster_estimates) {

            abundance_estimates_writer.addEstimates(path_cluster_estimates);
        }
    
        abundance_estimates_writer.close();
    }    

    double time_end = gbwt::readTimer();
    cerr << "Inferred path posterior probabilities" << ((inference_model != "haplotypes") ? " and abundances" : "") << " (" << time_end - time_clust << " seconds, " << gbwt::inGigabytes(gbwt::memoryUsage()) << " GB)" << endl;

	return 0;
}

