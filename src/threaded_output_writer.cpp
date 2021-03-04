
#include "threaded_output_writer.hpp"

#include <iomanip>

ThreadedOutputWriter::ThreadedOutputWriter(const string & filename, const string & compression_mode, const uint32_t num_threads) {

    writer_stream = bgzf_open(filename.c_str(), compression_mode.c_str());

    output_queue = new ProducerConsumerQueue<stringstream *>(num_threads * 5);
    writing_thread = thread(&ThreadedOutputWriter::write, this);
}

void ThreadedOutputWriter::close() {

    output_queue->pushedLast();

    writing_thread.join();
    delete output_queue;

    assert(bgzf_close(writer_stream) == 0);
}

void ThreadedOutputWriter::write() {

    stringstream * out_sstream = nullptr;

    while (output_queue->pop(&out_sstream)) {

        const string & tmp_out_string = out_sstream->str();   
        assert(bgzf_write(writer_stream, tmp_out_string.data(), tmp_out_string.size()) >= 0);

        delete out_sstream;
    }
}


ProbabilityClusterWriter::ProbabilityClusterWriter(const string filename_prefix, const uint32_t num_threads, const double prob_precision_in) : ThreadedOutputWriter(filename_prefix + "_probs.txt.gz", "wg", num_threads), prob_precision(prob_precision_in), prob_precision_digits(ceil(-1 * log10(prob_precision))) {}

void ProbabilityClusterWriter::addCluster(const vector<ReadPathProbabilities> & read_path_cluster_probs, const vector<PathInfo> & cluster_paths) {

    assert(!cluster_paths.empty());

    if (!read_path_cluster_probs.empty()) {

        auto out_sstream = new stringstream;

        *out_sstream << "#" << endl;
        *out_sstream << setprecision(3);
        *out_sstream << cluster_paths.front().name << "," << cluster_paths.front().length << "," << cluster_paths.front().effective_length;

        for (size_t i = 1; i < cluster_paths.size(); ++i) {

            *out_sstream << " " << cluster_paths.at(i).name << "," << cluster_paths.at(i).length << "," << cluster_paths.at(i).effective_length;
        }

        *out_sstream << endl;

        if (!read_path_cluster_probs.empty()) {

            *out_sstream << setprecision(prob_precision_digits);

            for (auto & read_path_probs: read_path_cluster_probs) {

                *out_sstream << read_path_probs.readCount() << " " << read_path_probs.noiseProb();

                for (auto & path_probs: read_path_probs.pathProbs()) {

                    *out_sstream << " " << path_probs.first << ":";

                    bool is_first = true;

                    for (auto path: path_probs.second) {

                        if (is_first) {

                            *out_sstream << path;
                            is_first = false;

                        } else {

                            *out_sstream << "," << path;
                        }
                    }
                }

                *out_sstream << endl;
            }
        }

        output_queue->push(out_sstream);
    }
}


GibbsSamplesWriter::GibbsSamplesWriter(const string filename_prefix, const uint32_t num_threads, const uint32_t num_gibbs_samples_in) : ThreadedOutputWriter(filename_prefix + "_gibbs.txt.gz", "wg", num_threads), num_gibbs_samples(num_gibbs_samples_in) {

    auto out_sstream = new stringstream;
    *out_sstream << "Name\tClusterID\tHaplotypeSampleId";

    for (size_t i = 1; i <= num_gibbs_samples; ++i) {

        *out_sstream << "\tReadCountSample_" << i;
    }

    *out_sstream << endl;
    output_queue->push(out_sstream);
}

void GibbsSamplesWriter::addSamples(const pair<uint32_t, PathClusterEstimates> & path_cluster_estimate) {

    if (!path_cluster_estimate.second.gibbs_read_count_samples.empty()) {

        uint32_t cur_hap_sample_id = 0;
        auto out_sstream = new stringstream;

        for (auto & read_count_samples: path_cluster_estimate.second.gibbs_read_count_samples) {

            assert(!read_count_samples.path_ids.empty());
            assert(read_count_samples.path_ids.size() == read_count_samples.samples.size());

            assert(read_count_samples.samples.front().size() % num_gibbs_samples == 0);

            for (size_t i = 0; i < read_count_samples.samples.front().size(); i += num_gibbs_samples) {

                ++cur_hap_sample_id;

                for (size_t j = 0; j < read_count_samples.path_ids.size(); ++j) {

                    assert(read_count_samples.samples.front().size() == read_count_samples.samples.at(j).size());

                    *out_sstream << path_cluster_estimate.second.paths.at(read_count_samples.path_ids.at(j)).name;
                    *out_sstream << "\t" << path_cluster_estimate.first;
                    *out_sstream << "\t" << cur_hap_sample_id;

                    for (size_t k = 0; k < num_gibbs_samples; ++k) {

                        *out_sstream << "\t" << read_count_samples.samples.at(j).at(i + k);
                    }

                    *out_sstream << endl;
                }
            }
        }

        output_queue->push(out_sstream);
    }
}


PosteriorEstimatesWriter::PosteriorEstimatesWriter(const string filename_prefix, const uint32_t num_threads, const uint32_t ploidy_in, const double min_posterior_in) : ThreadedOutputWriter(filename_prefix + ".txt", "wu", num_threads), ploidy(ploidy_in), min_posterior(min_posterior_in) {

    auto out_sstream = new stringstream;

    for (uint32_t i = 0; i < ploidy; ++i) {

        *out_sstream << "Name" << i + 1 << "\t";
    }

    *out_sstream << "ClusterID\tPosteriorProbability" << endl;
    output_queue->push(out_sstream);
}

void PosteriorEstimatesWriter::addEstimates(const vector<pair<uint32_t, PathClusterEstimates> > & path_cluster_estimates) {

    auto out_sstream = new stringstream;

    for (auto & cur_estimates: path_cluster_estimates) {

        assert(cur_estimates.second.path_group_sets.size() == cur_estimates.second.posteriors.cols());

        for (size_t i = 0; i < cur_estimates.second.path_group_sets.size(); ++i) {

            assert(cur_estimates.second.path_group_sets.at(i).size() == ploidy);

            if (cur_estimates.second.posteriors(0, i) >= min_posterior) {

                for (auto & path_idx: cur_estimates.second.path_group_sets.at(i)) {

                    *out_sstream << cur_estimates.second.paths.at(path_idx).name << "\t";
                }

                *out_sstream << cur_estimates.first;
                *out_sstream << "\t" << cur_estimates.second.posteriors(0, i);
                *out_sstream << endl;
            }
        }
    }

    output_queue->push(out_sstream);
}


AbundanceEstimatesWriter::AbundanceEstimatesWriter(const string filename_prefix, const uint32_t num_threads, const double total_transcript_count_in) : ThreadedOutputWriter(filename_prefix + ".txt", "wu", num_threads), total_transcript_count(total_transcript_count_in) {

    auto out_sstream = new stringstream;
    *out_sstream << "Name\tClusterID\tLength\tEffectiveLength\tHaplotypeProbability\tReadCount\tTPM" << endl;
    output_queue->push(out_sstream);
}

void AbundanceEstimatesWriter::addEstimates(const vector<pair<uint32_t, PathClusterEstimates> > & path_cluster_estimates) {

    auto out_sstream = new stringstream;

    for (auto & cur_estimates: path_cluster_estimates) {

        for (size_t i = 0; i < cur_estimates.second.paths.size(); ++i) {

            double transcript_count = 0;

            if (cur_estimates.second.paths.at(i).effective_length > 0) {

                transcript_count = cur_estimates.second.abundances(0, i) / cur_estimates.second.paths.at(i).effective_length;
            }

            *out_sstream << cur_estimates.second.paths.at(i).name;
            *out_sstream << "\t" << cur_estimates.first;
            *out_sstream << "\t" << cur_estimates.second.paths.at(i).length;
            *out_sstream << "\t" << cur_estimates.second.paths.at(i).effective_length;
            *out_sstream << "\t" << cur_estimates.second.posteriors(0, i);
            *out_sstream << "\t" << cur_estimates.second.abundances(0, i);
            *out_sstream << "\t" << transcript_count / total_transcript_count * pow(10, 6);
            *out_sstream << endl;
        }
    }
    
    output_queue->push(out_sstream);
}

