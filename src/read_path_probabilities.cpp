
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

double ReadPathProbabilities::noiseProb() const {

    return noise_prob;
}

const vector<pair<double, vector<uint32_t> > >& ReadPathProbabilities::pathProbs() const {

    return path_probs;
}

void ReadPathProbabilities::addReadCount(const uint32_t multiplicity_in) {

    read_count += multiplicity_in;
}

void ReadPathProbabilities::calcAlignPathProbs(const vector<AlignmentPath> & align_paths, const vector<vector<gbwt::size_type> > & align_paths_ids, const spp::sparse_hash_map<uint32_t, uint32_t> & clustered_path_index, const vector<PathInfo> & cluster_paths, const FragmentLengthDist & fragment_length_dist, const bool is_single_end, const double min_noise_prob) {

    assert(align_paths.size() > 1);
    assert(align_paths.size() == align_paths_ids.size());
    assert(clustered_path_index.size() == cluster_paths.size());

    assert(path_probs.empty());

    if (align_paths.front().min_mapq > 0) {

        noise_prob = max(prob_precision, max(min_noise_prob, Utils::phred_to_prob(align_paths.front().min_mapq)));
        assert(noise_prob < 1 && noise_prob > 0);

        assert(align_paths.back().gbwt_search.first.empty());
        assert(align_paths_ids.back().empty());
        assert(align_paths.back().score_sum <= 0);

        noise_prob += (1 - noise_prob) * exp(align_paths.back().score_sum * Utils::noise_score_log_base);

        if (align_paths.back().score_sum == 0) {

            assert(Utils::doubleCompare(noise_prob, 1));
            return;
        }

        vector<double> align_paths_log_probs;
        align_paths_log_probs.reserve(align_paths.size() - 1);

        for (size_t i = 0; i < align_paths.size() - 1; ++i) {

            const AlignmentPath & align_path = align_paths.at(i);

            assert(align_paths.front().min_mapq == align_path.min_mapq);
            align_paths_log_probs.emplace_back(align_path.score_sum * Utils::score_log_base);

            if (!is_single_end) {

                align_paths_log_probs.back() += fragment_length_dist.logProb(align_path.frag_length);
            }
        }
        
        assert(align_paths_ids.size() == align_paths_log_probs.size() + 1);
        vector<double> read_path_log_probs(clustered_path_index.size(), numeric_limits<double>::lowest());

        for (size_t i = 0; i < align_paths_ids.size() - 1; ++i) {

            assert(!align_paths_ids.at(i).empty());

            for (auto path_id: align_paths_ids.at(i)) {

                auto clustered_path_index_it = clustered_path_index.find(path_id);
                assert(clustered_path_index_it != clustered_path_index.end());

                uint32_t path_idx = clustered_path_index_it->second;

                if (Utils::doubleCompare(cluster_paths.at(path_idx).effective_length, 0)) {

                    assert(Utils::doubleCompare(read_path_log_probs.at(path_idx), numeric_limits<double>::lowest()));
                    read_path_log_probs.at(path_idx) = numeric_limits<double>::lowest();

                } else {

                    // account for rare cases when a mpmap alignment can have multiple alignments on the same path
                    read_path_log_probs.at(path_idx) = max(read_path_log_probs.at(path_idx), align_paths_log_probs.at(i) - log(cluster_paths.at(path_idx).effective_length));
                }
            }
        }

        double read_path_log_probs_sum = numeric_limits<double>::lowest();

        for (auto & log_prob: read_path_log_probs) {

            read_path_log_probs_sum = Utils::add_log(read_path_log_probs_sum, log_prob);
        }

        assert(read_path_log_probs_sum > numeric_limits<double>::lowest());

        for (size_t i = 0; i < read_path_log_probs.size(); ++i) {

            read_path_log_probs.at(i) = exp(read_path_log_probs.at(i) - read_path_log_probs_sum);
            read_path_log_probs.at(i) *= (1 - noise_prob);

            if (read_path_log_probs.at(i) >= prob_precision) {

                auto path_probs_it = path_probs.begin();

                while (path_probs_it != path_probs.end()) {

                    if (abs(path_probs_it->first - read_path_log_probs.at(i)) < prob_precision) {

                        path_probs_it->first = ((path_probs_it->first * path_probs_it->second.size() + read_path_log_probs.at(i)) / (path_probs_it->second.size() + 1));
                        path_probs_it->second.emplace_back(i);

                        break;
                    }
                    
                    ++path_probs_it;
                }

                if (path_probs_it == path_probs.end()) {

                    path_probs.emplace_back(read_path_log_probs.at(i), vector<uint32_t>({static_cast<uint32_t>(i)}));
                }
            }
        }

        assert(!path_probs.empty());
        sort(path_probs.begin(), path_probs.end());
    }
}

bool ReadPathProbabilities::quickMergeIdentical(const ReadPathProbabilities & probs_2) {

    if (pathProbs().size() != probs_2.pathProbs().size()) {

        return false;
    }

    if (abs(noise_prob - probs_2.noiseProb()) < prob_precision) {

        for (size_t i = 0; i < pathProbs().size(); ++i) {

            if (abs(pathProbs().at(i).first - probs_2.pathProbs().at(i).first) >= prob_precision) {

                return false;
            }

            if (pathProbs().at(i).second != probs_2.pathProbs().at(i).second) {

                return false;
            }  
        }

        addReadCount(probs_2.readCount());
        return true;
    } 

    return false;
}

bool operator==(const ReadPathProbabilities & lhs, const ReadPathProbabilities & rhs) { 

    if (lhs.readCount() == rhs.readCount() && Utils::doubleCompare(lhs.noiseProb(), rhs.noiseProb())) {

        if (lhs.pathProbs().size() == rhs.pathProbs().size()) {

            for (size_t i = 0; i < lhs.pathProbs().size(); ++i) {

                if (!Utils::doubleCompare(lhs.pathProbs().at(i).first, rhs.pathProbs().at(i).first)) {

                    return false;
                }

                if (lhs.pathProbs().at(i).second != rhs.pathProbs().at(i).second) {

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

    if (!Utils::doubleCompare(lhs.noiseProb(), rhs.noiseProb())) {

        return (lhs.noiseProb() < rhs.noiseProb());    
    } 

    if (lhs.pathProbs().size() != rhs.pathProbs().size()) {

        return (lhs.pathProbs().size() < rhs.pathProbs().size());
    }

    for (size_t i = 0; i < lhs.pathProbs().size(); ++i) {

        if (!Utils::doubleCompare(lhs.pathProbs().at(i).first, rhs.pathProbs().at(i).first)) {

            return (lhs.pathProbs().at(i).first < rhs.pathProbs().at(i).first);    
        }      

        if (lhs.pathProbs().at(i).second.size() != rhs.pathProbs().at(i).second.size()) {

            return (lhs.pathProbs().at(i).second.size() < rhs.pathProbs().at(i).second.size());    
        }  

        for (size_t j = 0; j < lhs.pathProbs().at(i).second.size(); ++j) {

            if (lhs.pathProbs().at(i).second.at(j) != rhs.pathProbs().at(i).second.at(j)) {

                return (lhs.pathProbs().at(i).second.at(j) < rhs.pathProbs().at(i).second.at(j));    
            } 
        }   
    }

    if (lhs.readCount() != rhs.readCount()) {

        return (lhs.readCount() < rhs.readCount());
    }

    return false;
}

ostream & operator<<(ostream & os, const ReadPathProbabilities & read_path_probs) {

    os << read_path_probs.readCount() << " | " << read_path_probs.noiseProb() << " |";

    for (auto & path_probs: read_path_probs.pathProbs()) {

        os << " " << path_probs.first << ": " << path_probs.second;
    }

    return os;
}

