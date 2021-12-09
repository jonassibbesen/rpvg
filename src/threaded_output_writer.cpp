
#include "threaded_output_writer.hpp"

#include <iomanip>

const uint32_t out_precision_digits = 8;

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


ProbabilityClusterWriter::ProbabilityClusterWriter(const string filename_prefix, const uint32_t num_threads, const double prob_precision_in) : ThreadedOutputWriter(filename_prefix + ".txt.gz", "wg", num_threads), prob_precision(prob_precision_in), prob_precision_digits(max(out_precision_digits, static_cast<uint32_t>(ceil(-1 * log10(prob_precision))))) {}

void ProbabilityClusterWriter::addCluster(const vector<ReadPathProbabilities> & read_path_cluster_probs, const vector<PathInfo> & cluster_paths) {

    assert(!cluster_paths.empty());

    if (!read_path_cluster_probs.empty()) {

        auto out_sstream = new stringstream;

        *out_sstream << "#" << endl;
        *out_sstream << setprecision(out_precision_digits);
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


ReadCountGibbsSamplesWriter::ReadCountGibbsSamplesWriter(const string filename_prefix, const uint32_t num_threads, const uint32_t num_gibbs_samples_in) : ThreadedOutputWriter(filename_prefix + ".txt.gz", "wg", num_threads), num_gibbs_samples(num_gibbs_samples_in) {

    noise_counts = vector<double>(num_gibbs_samples, 0);

    auto out_sstream = new stringstream;
    *out_sstream << "Name\tClusterID";

    for (uint32_t i = 0; i < num_gibbs_samples; ++i) {

        *out_sstream << "\tReadCountSample_" << i + 1;
    }

    *out_sstream << endl;
    output_queue->push(out_sstream);
}

void ReadCountGibbsSamplesWriter::addSamples(const pair<uint32_t, PathClusterEstimates> & path_cluster_estimate) {

    if (!path_cluster_estimate.second.gibbs_read_count_samples.empty()) {

        uint32_t noise_count_idx = 0;

        vector<vector<uint32_t> > path_gibbs_sampling_index(path_cluster_estimate.second.paths.size());

        for (size_t i = 0; i < path_cluster_estimate.second.gibbs_read_count_samples.size(); ++i) {

            const CountSamples & count_samples = path_cluster_estimate.second.gibbs_read_count_samples.at(i);

            assert(!count_samples.path_ids.empty());
            assert(!count_samples.abundance_samples.empty());

            assert(count_samples.abundance_samples.size() % count_samples.path_ids.size() == 0);
            assert(count_samples.abundance_samples.size() / count_samples.path_ids.size() == count_samples.noise_samples.size());

            for (auto & noise_sample: count_samples.noise_samples) {

                assert(noise_sample < path_cluster_estimate.second.total_count || Utils::doubleCompare(noise_sample, path_cluster_estimate.second.total_count));

                noise_counts.at(noise_count_idx) += noise_sample;
                ++noise_count_idx;  
            }

            for (size_t j = 0; j < count_samples.path_ids.size(); ++j) {

                vector<uint32_t> * sampling_indices = &(path_gibbs_sampling_index.at(count_samples.path_ids.at(j)));

                if (sampling_indices->empty()) {

                    *sampling_indices = vector<uint32_t>(path_cluster_estimate.second.gibbs_read_count_samples.size(), numeric_limits<uint32_t>::max());
                }

                sampling_indices->at(i) = j;
            }            
        }

        while (noise_count_idx < num_gibbs_samples) {

            noise_counts.at(noise_count_idx) += path_cluster_estimate.second.total_count;
            ++noise_count_idx;                    
        }

        auto out_sstream = new stringstream;
        *out_sstream << setprecision(out_precision_digits);

        for (size_t i = 0; i < path_gibbs_sampling_index.size(); ++i) {

            const vector<uint32_t> & sampling_indices = path_gibbs_sampling_index.at(i);

            if (!sampling_indices.empty()) {

                *out_sstream << path_cluster_estimate.second.paths.at(i).name;
                *out_sstream << "\t" << path_cluster_estimate.first;

                uint32_t num_samples = 0;

                for (size_t j = 0; j < sampling_indices.size(); ++j) {

                    const CountSamples & count_samples = path_cluster_estimate.second.gibbs_read_count_samples.at(j);

                    if (sampling_indices.at(j) == numeric_limits<uint32_t>::max()) {

                        for (size_t k = 0; k < count_samples.abundance_samples.size() / count_samples.path_ids.size(); ++k) {

                            *out_sstream << "\t0";
                            ++num_samples;
                        }

                    } else {

                        for (size_t k = 0; k < count_samples.abundance_samples.size() / count_samples.path_ids.size(); ++k) {

                            *out_sstream << "\t" << count_samples.abundance_samples.at(k * count_samples.path_ids.size() + sampling_indices.at(j));
                            ++num_samples;
                        }
                    }
                }

                while (num_samples < num_gibbs_samples) {
    
                    *out_sstream << "\t0";
                    ++num_samples;                    
                }

                *out_sstream << endl;                
            }
        }

        output_queue->push(out_sstream);
    
    } else {

        for (auto & noise_count: noise_counts) {

            noise_count += path_cluster_estimate.second.total_count;
        }
    }
}

void ReadCountGibbsSamplesWriter::addNoiseTranscript(const uint32_t unaligned_read_count) {
    
    auto out_sstream = new stringstream;
    *out_sstream << setprecision(out_precision_digits);

    *out_sstream << "Unknown\t0";

    for (auto & noise_count: noise_counts) {

        *out_sstream << "\t" << noise_count + unaligned_read_count;
    }

    *out_sstream << endl;
    output_queue->push(out_sstream);
}


JointHaplotypeEstimatesWriter::JointHaplotypeEstimatesWriter(const string filename_prefix, const uint32_t num_threads, const uint32_t ploidy_in, const double min_posterior_in) : ThreadedOutputWriter(filename_prefix + ".txt", "wu", num_threads), ploidy(ploidy_in), min_posterior(min_posterior_in) {

    auto out_sstream = new stringstream;

    for (uint32_t i = 0; i < ploidy; ++i) {

        *out_sstream << "Name_" << i + 1 << "\t";
    }

    *out_sstream << "ClusterID\tHaplotypingProbability" << endl;
    output_queue->push(out_sstream);
}

void JointHaplotypeEstimatesWriter::addEstimates(const vector<pair<uint32_t, PathClusterEstimates> > & path_cluster_estimates) {

    auto out_sstream = new stringstream;
    *out_sstream << setprecision(out_precision_digits);

    for (auto & cur_estimates: path_cluster_estimates) {

        assert(cur_estimates.second.posteriors.size() == cur_estimates.second.path_group_sets.size());

        for (size_t i = 0; i < cur_estimates.second.path_group_sets.size(); ++i) {

            assert(!cur_estimates.second.path_group_sets.at(i).empty());
            assert(cur_estimates.second.path_group_sets.at(i).size() <= ploidy);

            if (cur_estimates.second.posteriors.at(i) >= min_posterior) {

                for (auto & path: cur_estimates.second.path_group_sets.at(i)) {

                    *out_sstream << cur_estimates.second.paths.at(path).name << "\t";
                }

                for (size_t j = cur_estimates.second.path_group_sets.at(i).size(); j < ploidy; ++j) {

                    *out_sstream << ".\t";
                }

                *out_sstream << cur_estimates.first;
                *out_sstream << "\t" << cur_estimates.second.posteriors.at(i);
                *out_sstream << endl;
            }
        }
    }

    output_queue->push(out_sstream);
}


AbundanceEstimatesWriter::AbundanceEstimatesWriter(const string filename_prefix, const uint32_t num_threads, const double total_transcript_count_in) : ThreadedOutputWriter(filename_prefix + ".txt", "wu", num_threads), total_transcript_count(total_transcript_count_in) {

    noise_count = 0;

    auto out_sstream = new stringstream;
    *out_sstream << setprecision(out_precision_digits);

    *out_sstream << "Name\tClusterID\tLength\tEffectiveLength\tReadCount\tTPM" << endl;
    output_queue->push(out_sstream);
}

void AbundanceEstimatesWriter::addEstimates(const vector<pair<uint32_t, PathClusterEstimates> > & path_cluster_estimates) {

    auto out_sstream = new stringstream;
    *out_sstream << setprecision(out_precision_digits);

    for (auto & cur_estimates: path_cluster_estimates) {

        assert(cur_estimates.second.paths.size() == cur_estimates.second.path_group_sets.size());
        assert(cur_estimates.second.paths.size() == cur_estimates.second.abundances.size());

        double abundance_sum = 0;

        for (size_t i = 0; i < cur_estimates.second.path_group_sets.size(); ++i) {

            assert(cur_estimates.second.path_group_sets.at(i).size() == 1);
            auto path = cur_estimates.second.path_group_sets.at(i).front();

            double transcript_count = 0;

            if (cur_estimates.second.paths.at(path).effective_length > 0) {

                transcript_count = cur_estimates.second.abundances.at(i) / cur_estimates.second.paths.at(path).effective_length;
            }

            *out_sstream << cur_estimates.second.paths.at(path).name;
            *out_sstream << "\t" << cur_estimates.first;
            *out_sstream << "\t" << cur_estimates.second.paths.at(path).length;
            *out_sstream << "\t" << cur_estimates.second.paths.at(path).effective_length;
            *out_sstream << "\t" << cur_estimates.second.abundances.at(i);
            *out_sstream << "\t" << transcript_count / total_transcript_count * pow(10, 6);
            *out_sstream << endl;

            abundance_sum += cur_estimates.second.abundances.at(i);
        }

        assert(cur_estimates.second.noise_count < cur_estimates.second.total_count || Utils::doubleCompare(cur_estimates.second.noise_count, cur_estimates.second.total_count));
        assert(Utils::doubleCompare(abundance_sum + cur_estimates.second.noise_count, cur_estimates.second.total_count));

        noise_count += cur_estimates.second.noise_count;
    }
    
    output_queue->push(out_sstream);
}

void AbundanceEstimatesWriter::addNoiseTranscript(const uint32_t unaligned_read_count) {
    
    auto out_sstream = new stringstream;
    *out_sstream << setprecision(out_precision_digits);

    *out_sstream << "Unknown\t0\t0\t0\t" << noise_count + unaligned_read_count << "\t0" << endl;
    output_queue->push(out_sstream);
}


HaplotypeAbundanceEstimatesWriter::HaplotypeAbundanceEstimatesWriter(const string filename_prefix, const uint32_t num_threads, const uint32_t ploidy_in, const double total_transcript_count_in) : ThreadedOutputWriter(filename_prefix + ".txt", "wu", num_threads), ploidy(ploidy_in), total_transcript_count(total_transcript_count_in) {

    noise_count = 0;

    auto out_sstream = new stringstream;
    *out_sstream << "Name\tClusterID\tLength\tEffectiveLength\tHaplotypeProbability\tReadCount\tTPM" << endl;
    output_queue->push(out_sstream);
}

void HaplotypeAbundanceEstimatesWriter::addEstimates(const vector<pair<uint32_t, PathClusterEstimates> > & path_cluster_estimates) {

    auto out_sstream = new stringstream;
    *out_sstream << setprecision(out_precision_digits);

    for (auto & cur_estimates: path_cluster_estimates) {

        assert(cur_estimates.second.path_group_sets.size() == cur_estimates.second.posteriors.size());

        vector<double> haplotype_probs(cur_estimates.second.paths.size(), 0);
        vector<double> read_counts(cur_estimates.second.paths.size(), 0);

        auto abundances_it = cur_estimates.second.abundances.begin();

        for (size_t i = 0; i < cur_estimates.second.path_group_sets.size(); ++i) {

            assert(!cur_estimates.second.path_group_sets.at(i).empty());
            assert(cur_estimates.second.path_group_sets.at(i).size() <= ploidy);

            haplotype_probs.at(cur_estimates.second.path_group_sets.at(i).front()) += cur_estimates.second.posteriors.at(i);
            
            read_counts.at(cur_estimates.second.path_group_sets.at(i).front()) += *abundances_it;
            ++abundances_it;

            for (size_t j = 1; j < cur_estimates.second.path_group_sets.at(i).size(); ++j) {

                if (cur_estimates.second.path_group_sets.at(i).at(j) != cur_estimates.second.path_group_sets.at(i).at(j - 1)) {

                    haplotype_probs.at(cur_estimates.second.path_group_sets.at(i).at(j)) += cur_estimates.second.posteriors.at(i);
                }

                read_counts.at(cur_estimates.second.path_group_sets.at(i).at(j)) += *abundances_it;
                ++abundances_it;
            }
        }

        assert(abundances_it == cur_estimates.second.abundances.end());

        double abundance_sum = 0;

        for (size_t i = 0; i < cur_estimates.second.paths.size(); ++i) {

            double transcript_count = 0;

            if (cur_estimates.second.paths.at(i).effective_length > 0) {

                transcript_count = read_counts.at(i) / cur_estimates.second.paths.at(i).effective_length;
            }

            *out_sstream << cur_estimates.second.paths.at(i).name;
            *out_sstream << "\t" << cur_estimates.first;
            *out_sstream << "\t" << cur_estimates.second.paths.at(i).length;
            *out_sstream << "\t" << cur_estimates.second.paths.at(i).effective_length;
            *out_sstream << "\t" << haplotype_probs.at(i);
            *out_sstream << "\t" << read_counts.at(i);
            *out_sstream << "\t" << transcript_count / total_transcript_count * pow(10, 6);
            *out_sstream << endl;

            abundance_sum += read_counts.at(i);
        }

        assert(cur_estimates.second.noise_count < cur_estimates.second.total_count || Utils::doubleCompare(cur_estimates.second.noise_count, cur_estimates.second.total_count));
        assert(Utils::doubleCompare(abundance_sum + cur_estimates.second.noise_count, cur_estimates.second.total_count));

        noise_count += cur_estimates.second.noise_count;
    }
    
    output_queue->push(out_sstream);
}

void HaplotypeAbundanceEstimatesWriter::addNoiseTranscript(const uint32_t unaligned_read_count) {
    
    auto out_sstream = new stringstream;
    *out_sstream << setprecision(out_precision_digits);

    *out_sstream << "Unknown\t0\t0\t0\t0\t" << noise_count + unaligned_read_count << "\t0" << endl;
    output_queue->push(out_sstream);
}

JointHaplotypeAbundanceEstimatesWriter::JointHaplotypeAbundanceEstimatesWriter(const string filename_prefix, const uint32_t num_threads, const uint32_t ploidy_in, const double min_posterior_in, const double total_transcript_count_in) : ThreadedOutputWriter(filename_prefix + ".txt", "wu", num_threads), ploidy(ploidy_in), min_posterior(min_posterior_in), total_transcript_count(total_transcript_count_in) {

    noise_counts = vector<double>(ploidy, 0);

    auto out_sstream = new stringstream;

    for (uint32_t i = 0; i < ploidy; ++i) {

        *out_sstream << "Name_" << i + 1 << "\t";
    }

    *out_sstream << "ClusterID\tHaplotypingProbability";

    for (uint32_t i = 0; i < ploidy; ++i) {

        *out_sstream << "\tReadCount_" << i + 1 << "\tTPM_" << i + 1;
    }

    *out_sstream << endl;
    output_queue->push(out_sstream);
}

void JointHaplotypeAbundanceEstimatesWriter::addEstimates(const vector<pair<uint32_t, PathClusterEstimates> > & path_cluster_estimates) {

    auto out_sstream = new stringstream;
    *out_sstream << setprecision(out_precision_digits);

    for (auto & cur_estimates: path_cluster_estimates) {

        assert(cur_estimates.second.posteriors.size() == cur_estimates.second.path_group_sets.size());

        auto abundances_it = cur_estimates.second.abundances.begin();
        double abundance_sum = 0;

        for (size_t i = 0; i < cur_estimates.second.path_group_sets.size(); ++i) {

            assert(!cur_estimates.second.path_group_sets.at(i).empty());
            assert(cur_estimates.second.path_group_sets.at(i).size() <= ploidy);

            if (cur_estimates.second.posteriors.at(i) >= min_posterior) {

                for (auto & path: cur_estimates.second.path_group_sets.at(i)) {

                    *out_sstream << cur_estimates.second.paths.at(path).name << "\t";
                }

                for (size_t j = cur_estimates.second.path_group_sets.at(i).size(); j < ploidy; ++j) {

                    *out_sstream << ".\t";
                }

                *out_sstream << cur_estimates.first;
                *out_sstream << "\t" << cur_estimates.second.posteriors.at(i);

                for (auto & path: cur_estimates.second.path_group_sets.at(i)) {

                    *out_sstream << "\t" << *abundances_it;
                    abundance_sum += *abundances_it;

                    double transcript_count = 0;

                    if (cur_estimates.second.paths.at(path).effective_length > 0) {

                        transcript_count = *abundances_it / cur_estimates.second.paths.at(path).effective_length;
                    }

                    *out_sstream << "\t" << transcript_count / total_transcript_count * pow(10, 6);
                    ++abundances_it;
                }

                for (size_t j = cur_estimates.second.path_group_sets.at(i).size(); j < ploidy; ++j) {

                    *out_sstream << "\t0\t0";
                }

                *out_sstream << endl;
            }
        }

        assert(abundances_it == cur_estimates.second.abundances.end());

        assert(cur_estimates.second.noise_count < cur_estimates.second.total_count || Utils::doubleCompare(cur_estimates.second.noise_count, cur_estimates.second.total_count));
        assert(Utils::doubleCompare(abundance_sum + cur_estimates.second.noise_count, cur_estimates.second.total_count));

        for (auto & noise_count: noise_counts) {

            noise_count += cur_estimates.second.noise_count / noise_counts.size();
        }
    }

    output_queue->push(out_sstream);
}

void JointHaplotypeAbundanceEstimatesWriter::addNoiseTranscript(const uint32_t unaligned_read_count) {
    
    auto out_sstream = new stringstream;
    *out_sstream << setprecision(out_precision_digits);

    for (uint32_t i = 0; i < ploidy; ++i) {

        *out_sstream << "Unknown\t";
    }

    *out_sstream << "0\t0";

    for (auto & noise_count: noise_counts) {

        *out_sstream << "\t" << noise_count + static_cast<float>(unaligned_read_count) / noise_counts.size() << "\t0";
    }

    *out_sstream << endl;
    output_queue->push(out_sstream);
}

