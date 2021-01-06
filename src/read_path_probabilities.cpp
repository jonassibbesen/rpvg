
#include "read_path_probabilities.hpp"

#include <assert.h>
#include <algorithm>
#include <numeric>
#include <limits>
#include <sstream>


ReadPathProbabilities::ReadPathProbabilities() {

    read_count = 0;
    noise_prob = 1;

    prob_precision = pow(10, -8);
}

ReadPathProbabilities::ReadPathProbabilities(const uint32_t read_count_in, const double prob_precision_in) : read_count(read_count_in), prob_precision(prob_precision_in) {
    
    noise_prob = 1;
}

uint32_t ReadPathProbabilities::readCount() const {

    return read_count;
}

double ReadPathProbabilities::noiseProbability() const {

    return noise_prob;
}

const vector<pair<uint32_t, double> > & ReadPathProbabilities::probabilities() const {

    return read_path_probs;
}

void ReadPathProbabilities::addReadCount(const uint32_t multiplicity_in) {

    read_count += multiplicity_in;
}

void ReadPathProbabilities::calcReadPathProbabilities(const vector<AlignmentPath> & align_paths, const vector<vector<gbwt::size_type> > & align_paths_ids, const spp::sparse_hash_map<uint32_t, uint32_t> & clustered_path_index, const vector<PathInfo> & cluster_paths, const FragmentLengthDist & fragment_length_dist, const bool is_single_end) {

    assert(!align_paths.empty());
    assert(align_paths.size() == align_paths_ids.size());

    assert(clustered_path_index.size() == cluster_paths.size());

    if (align_paths.front().min_mapq > 0) {

        noise_prob = phred_to_prob(align_paths.front().min_mapq);
        assert(noise_prob < 1 && noise_prob > 0);

        vector<double> align_paths_log_probs;
        align_paths_log_probs.reserve(align_paths.size());

        for (auto & align_path: align_paths) {

            assert(align_paths.front().min_mapq == align_path.min_mapq);
            align_paths_log_probs.emplace_back(score_log_base * align_path.score_sum);

            if (!is_single_end) {

                align_paths_log_probs.back() += fragment_length_dist.logProb(align_path.frag_length);
            }
        }
        
        vector<double> read_path_log_probs(clustered_path_index.size(), numeric_limits<double>::lowest());

        for (size_t i = 0; i < align_paths.size(); ++i) {

            for (auto path_id: align_paths_ids.at(i)) {

                auto clustered_path_index_it = clustered_path_index.find(path_id);
                assert(clustered_path_index_it != clustered_path_index.end());

                uint32_t path_idx = clustered_path_index_it->second;

                if (doubleCompare(cluster_paths.at(path_idx).effective_length, 0)) {

                    assert(doubleCompare(read_path_log_probs.at(path_idx), numeric_limits<double>::lowest()));
                    read_path_log_probs.at(path_idx) = numeric_limits<double>::lowest();

                } else {

                    // account for really rare cases when a mpmap alignment can have multiple alignments on the same path
                    read_path_log_probs.at(path_idx) = max(read_path_log_probs.at(path_idx), align_paths_log_probs.at(i) - log(cluster_paths.at(path_idx).effective_length));
                }
            }
        }

        double read_path_log_probs_sum = numeric_limits<double>::lowest();

        for (auto & log_prob: read_path_log_probs) {

            read_path_log_probs_sum = add_log(read_path_log_probs_sum, log_prob);
        }

        assert(read_path_log_probs_sum > numeric_limits<double>::lowest());

        for (size_t i = 0; i < read_path_log_probs.size(); ++i) {

            read_path_log_probs.at(i) = exp(read_path_log_probs.at(i) - read_path_log_probs_sum);
            read_path_log_probs.at(i) *= (1 - noise_prob);

            if (read_path_log_probs.at(i) >= prob_precision) {

                read_path_probs.emplace_back(i, read_path_log_probs.at(i));
            }
        }
    }

    sort(read_path_probs.begin(), read_path_probs.end());
}

bool ReadPathProbabilities::mergeIdenticalReadPathProbabilities(const ReadPathProbabilities & probs_2) {

    if (probabilities().size() != probs_2.probabilities().size()) {

        return false;
    }

    if (abs(noise_prob - probs_2.noiseProbability()) < prob_precision) {

        for (size_t i = 0; i < read_path_probs.size(); ++i) {

            if (read_path_probs.at(i).first != probs_2.probabilities().at(i).first) {

                return false;
            }  

            if (abs(read_path_probs.at(i).second - probs_2.probabilities().at(i).second) >= prob_precision) {

                return false;
            }
        }

        addReadCount(probs_2.readCount());
        return true;
    } 

    return false;
}

vector<pair<double, vector<uint32_t> > > ReadPathProbabilities::collapsedProbabilities() const {

    vector<pair<double, vector<uint32_t> > > collapsed_probs;

    for (auto & prob: read_path_probs) {

        auto collapsed_probs_it = collapsed_probs.begin();

        while (collapsed_probs_it != collapsed_probs.end()) {

            if (abs(collapsed_probs_it->first - prob.second) < prob_precision) {

                collapsed_probs_it->second.emplace_back(prob.first);
                break;
            }
            
            ++collapsed_probs_it;
        }

        if (collapsed_probs_it == collapsed_probs.end()) {

            collapsed_probs.emplace_back(prob.second, vector<uint32_t>(1, prob.first));
        }
    }

    sort(collapsed_probs.begin(), collapsed_probs.end());
    return collapsed_probs;
}

bool operator==(const ReadPathProbabilities & lhs, const ReadPathProbabilities & rhs) { 

    if (lhs.readCount() == rhs.readCount() && doubleCompare(lhs.noiseProbability(), rhs.noiseProbability())) {

        if (lhs.probabilities().size() == rhs.probabilities().size()) {

            for (size_t i = 0; i < lhs.probabilities().size(); ++i) {

                if (lhs.probabilities().at(i).first != rhs.probabilities().at(i).first) {

                    return false;
                }  

                if (!doubleCompare(lhs.probabilities().at(i).second, rhs.probabilities().at(i).second)) {

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

    if (lhs.probabilities().size() != rhs.probabilities().size()) {

        return (lhs.probabilities().size() < rhs.probabilities().size());
    }

    for (size_t i = 0; i < lhs.probabilities().size(); ++i) {

        if (lhs.probabilities().at(i).first != rhs.probabilities().at(i).first) {

            return (lhs.probabilities().at(i).first < rhs.probabilities().at(i).first);    
        }  

        if (!doubleCompare(lhs.probabilities().at(i).second, rhs.probabilities().at(i).second)) {

            return (lhs.probabilities().at(i).second < rhs.probabilities().at(i).second);    
        }         
    }

    if (lhs.readCount() != rhs.readCount()) {

        return (lhs.readCount() < rhs.readCount());
    }

    return false;
}

ostream & operator<<(ostream & os, const ReadPathProbabilities & read_path_probs) {

    os << read_path_probs.readCount() << " | " << read_path_probs.noiseProbability() << " |";

    for (auto & prob: read_path_probs.probabilities()) {

        os << " " << prob.first << "," << prob.second;
    }

    return os;
}

