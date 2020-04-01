
#include "read_path_probabilities.hpp"

#include <assert.h>
#include <algorithm>
#include <numeric>
#include <limits>
#include <sstream>


ReadPathProbabilities::ReadPathProbabilities() : score_log_base(1) {

    noise_prob = 1;     
}

ReadPathProbabilities::ReadPathProbabilities(const uint32_t num_paths, const double score_log_base_in) : score_log_base(score_log_base_in) {

    noise_prob = 1;
    read_path_probs = vector<double>(num_paths, 0);
}

void ReadPathProbabilities::calcReadPathProbabilities(const vector<AlignmentPath> & align_paths, const unordered_map<uint32_t, uint32_t> & clustered_path_index, const FragmentLengthDist & fragment_length_dist, const bool is_single_end) {

    assert(!align_paths.empty());
    assert(clustered_path_index.size() == read_path_probs.size());

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

                read_path_probs.at(clustered_path_index.at(path)) = exp(align_paths_log_probs.at(i));
            }

            read_path_probs_sum += exp(align_paths_log_probs.at(i)) * align_paths.at(i).ids.size();
        }

        assert(read_path_probs_sum > 0);

        for (auto & prob: read_path_probs) {

            prob /= read_path_probs_sum;
            prob *= (1 - noise_prob);
        }
    }
}

void ReadPathProbabilities::addPositionalProbabilities(const vector<double> & path_lengths) {

    assert(path_lengths.size() == read_path_probs.size());

    if (noise_prob < 1) {

        double read_path_probs_sum = 0;

        for (size_t i = 0; i < read_path_probs.size(); ++i) {

            if (doubleCompare(path_lengths.at(i), 0)) {

                read_path_probs.at(i) = 0;

            } else {

                read_path_probs.at(i) /= path_lengths.at(i);
            }

            read_path_probs_sum += read_path_probs.at(i);
        }

        assert(read_path_probs_sum > 0);

        for (auto & probs: read_path_probs) {

            probs /= read_path_probs_sum;
            probs *= (1 - noise_prob);
        }  
    }
}

string ReadPathProbabilities::getCollapsedProbabilityString(const double precision) const {

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

    stringstream collpased_probs_ss;

    collpased_probs_ss << noise_prob;

    for (auto & prob: collpased_probs) {

        collpased_probs_ss << " " << prob.first << ":";
        bool is_first = true;

        for (auto idx: prob.second) {

            if (is_first) {

                collpased_probs_ss << idx;
                is_first = false;

            } else {

                collpased_probs_ss << "," << idx;
            }
        }
    }

    return collpased_probs_ss.str();
}

bool operator==(const ReadPathProbabilities & lhs, const ReadPathProbabilities & rhs) { 

    if (doubleCompare(lhs.noise_prob, rhs.noise_prob)) {

        if (lhs.read_path_probs.size() == rhs.read_path_probs.size()) {

            for (size_t i = 0; i < lhs.read_path_probs.size(); ++i) {

                if (!doubleCompare(lhs.read_path_probs.at(i), rhs.read_path_probs.at(i))) {

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

    if (!doubleCompare(lhs.noise_prob, rhs.noise_prob)) {

        return (lhs.noise_prob < rhs.noise_prob);    
    } 

    assert(lhs.read_path_probs.size() == rhs.read_path_probs.size());

    for (size_t i = 0; i < lhs.read_path_probs.size(); ++i) {

        if (!doubleCompare(lhs.read_path_probs.at(i), rhs.read_path_probs.at(i))) {

            return (lhs.read_path_probs.at(i) < rhs.read_path_probs.at(i));    
        }         
    }   

    return false;
}

ostream & operator<<(ostream & os, const ReadPathProbabilities & read_path_probs) {

    os << read_path_probs.noise_prob;

    for (auto & prob: read_path_probs.read_path_probs) {

        os << " " << prob;
    }

    return os;
}

