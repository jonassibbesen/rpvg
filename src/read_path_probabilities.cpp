
#include "read_path_probabilities.hpp"

#include <assert.h>
#include <algorithm>
#include <numeric>
#include <limits>
#include <sstream>


ReadPathProbabilities::ReadPathProbabilities() : score_log_base(1), fragment_length_dist(FragmentLengthDist()) {

    noise_prob = 1;     
}

ReadPathProbabilities::ReadPathProbabilities(const uint32_t num_paths, const double score_log_base_in, const FragmentLengthDist & fragment_length_dist_in) : score_log_base(score_log_base_in), fragment_length_dist(fragment_length_dist_in) {

    noise_prob = 1;
    read_path_probs = vector<double>(num_paths, 0);
}

const vector<double> & ReadPathProbabilities::probabilities() const {

    return read_path_probs;
}

double ReadPathProbabilities::noiseProbability() const {

    return noise_prob;
}

void ReadPathProbabilities::calcReadPathProbabilities(const vector<AlignmentPath> & align_paths, const unordered_map<uint32_t, uint32_t> & clustered_path_index, const vector<Path> & cluster_paths, const bool is_single_end) {

    assert(!align_paths.empty());
    assert(clustered_path_index.size() == read_path_probs.size());
    assert(cluster_paths.size() == read_path_probs.size());

    if (align_paths.front().mapq_comb > 0) {

        noise_prob = phred_to_prob(align_paths.front().mapq_comb);
        assert(noise_prob < 1);

        vector<double> align_paths_log_probs;
        align_paths_log_probs.reserve(align_paths.size());

        double align_paths_log_probs_sum = numeric_limits<double>::lowest();

        for (auto & align_path: align_paths) {

            align_paths_log_probs.emplace_back(score_log_base * align_path.score_sum);

            if (!is_single_end) {

                align_paths_log_probs.back() += fragment_length_dist.logProb(align_path.seq_length);
            }

            align_paths_log_probs_sum = add_log(align_paths_log_probs_sum, align_paths_log_probs.back());
        }

        for (auto & log_probs: align_paths_log_probs) {

            log_probs -= align_paths_log_probs_sum;
        }

        double read_path_probs_sum = 0;

        for (size_t i = 0; i < align_paths.size(); ++i) {

            for (auto & path: align_paths.at(i).ids) {

                uint32_t path_idx = clustered_path_index.at(path);

                read_path_probs.at(path_idx) = exp(align_paths_log_probs.at(i));
            
                if (doubleCompare(cluster_paths.at(path_idx).effective_length, 0)) {

                    read_path_probs.at(path_idx) = 0;

                } else {

                    read_path_probs.at(path_idx) /= cluster_paths.at(path_idx).effective_length;
                }

                read_path_probs_sum += read_path_probs.at(clustered_path_index.at(path));
            }
        }

        assert(read_path_probs_sum > 0);

        for (auto & prob: read_path_probs) {

            prob /= read_path_probs_sum;
            prob *= (1 - noise_prob);
        }
    }
}

vector<pair<double, vector<uint32_t> > > ReadPathProbabilities::collapsedProbabilities(const double precision) const {

    vector<pair<double, vector<uint32_t> > > collpased_probs;

    for (size_t i = 0; i < read_path_probs.size(); ++i) {

        auto collpased_probs_it = collpased_probs.begin();

        while (collpased_probs_it != collpased_probs.end()) {

            if (abs(collpased_probs_it->first - read_path_probs.at(i)) < precision) {

                collpased_probs_it->second.emplace_back(i);
                break;
            }
            
            ++collpased_probs_it;
        }

        if (collpased_probs_it == collpased_probs.end()) {

            collpased_probs.emplace_back(read_path_probs.at(i), vector<uint32_t>(1, i));
        }
    }

    return collpased_probs;
}

bool operator==(const ReadPathProbabilities & lhs, const ReadPathProbabilities & rhs) { 

    if (doubleCompare(lhs.noiseProbability(), rhs.noiseProbability())) {

        if (lhs.probabilities().size() == rhs.probabilities().size()) {

            for (size_t i = 0; i < lhs.probabilities().size(); ++i) {

                if (!doubleCompare(lhs.probabilities().at(i), rhs.probabilities().at(i))) {

                    return false;
                }
            }

            return true;
        }
    } 

    return false;
}

bool operator!=(const ReadPathProbabilities & lhs, const ReadPathProbabilities & rhs) { 

    return !(lhs == rhs);
}

bool operator<(const ReadPathProbabilities & lhs, const ReadPathProbabilities & rhs) { 

    if (!doubleCompare(lhs.noiseProbability(), rhs.noiseProbability())) {

        return (lhs.noiseProbability() < rhs.noiseProbability());    
    } 

    assert(lhs.probabilities().size() == rhs.probabilities().size());

    for (size_t i = 0; i < lhs.probabilities().size(); ++i) {

        if (!doubleCompare(lhs.probabilities().at(i), rhs.probabilities().at(i))) {

            return (lhs.probabilities().at(i) < rhs.probabilities().at(i));    
        }         
    }   

    return false;
}

ostream & operator<<(ostream & os, const ReadPathProbabilities & read_path_probs) {

    os << read_path_probs.noiseProbability();

    for (auto & prob: read_path_probs.probabilities()) {

        os << " " << prob;
    }

    return os;
}

