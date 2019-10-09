
#include "read_path_probs.hpp"

#include <assert.h>
#include <algorithm>
#include <numeric>

#include "gssw.h"


ReadPathProbs::ReadPathProbs() {

    noise_prob = 1;     
    score_log_base = 1;
}

ReadPathProbs::ReadPathProbs(const int32_t num_paths) {

    noise_prob = 1;
    read_path_probs = vector<double>(num_paths, 0);
    score_log_base = gssw_dna_recover_log_base(1, 4, 0.5, double_precision);
}

void ReadPathProbs::calcReadPathProbs(const vector<AlignmentPath> & align_paths, const unordered_map<int32_t, int32_t> & clustered_path_index, const FragmentLengthDist & fragment_length_dist) {

    assert(!align_paths.empty());
    assert(clustered_path_index.size() == read_path_probs.size());

    if (align_paths.front().mapqMin() > 0) {

        assert(align_paths.front().mapqs.size() == 2);
        assert(align_paths.front().scores.size() == 2);

        noise_prob = align_paths.front().mapqProb();
        assert(noise_prob < 1);

        auto align_paths_log_probs = calcRelativeAlignmentScoreLogProbs(align_paths);
        double align_paths_log_probs_sum = numeric_limits<double>::lowest();

        for (size_t i = 0; i < align_paths.size(); ++i) {

            assert(align_paths.at(i).mapqs.size() == 2);
            assert(align_paths.at(i).scores.size() == 2);

            align_paths_log_probs.at(i) += fragment_length_dist.logProb(align_paths.at(i).seq_length);
            align_paths_log_probs_sum = add_log(align_paths_log_probs_sum, align_paths_log_probs.at(i));
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

        for (auto & probs: read_path_probs) {

            probs /= read_path_probs_sum;
        }
    }
}

vector<double> ReadPathProbs::calcRelativeAlignmentScoreLogProbs(const vector<AlignmentPath> & align_paths) const {

    vector<double> align_score_log_probs;
    align_score_log_probs.reserve(align_paths.size());

    double score_log_probs_sum = numeric_limits<double>::lowest();

    for (auto & align_path: align_paths) {

        align_score_log_probs.emplace_back(score_log_base * align_path.scoreSum());
        score_log_probs_sum = add_log(score_log_probs_sum, align_score_log_probs.back());
    }

    for (auto & log_probs: align_score_log_probs) {

        log_probs -= score_log_probs_sum;
    }

    return align_score_log_probs;
}

double ReadPathProbs::calcReadMappingProbs(const vg::Alignment & alignment, const vector<double> & quality_match_probs, const vector<double> & quality_mismatch_probs, const double indel_prob) const {

    double align_path_prob = 0;

    auto & base_qualities = alignment.quality();
    int32_t cur_pos = 0;

    for (auto & mapping: alignment.path().mapping()) {

        for (auto & edit: mapping.edit()) {

            if (edit.from_length() == edit.to_length() && edit.sequence().empty()) {

                for (int32_t i = cur_pos; i < cur_pos + edit.from_length(); ++i) {

                    align_path_prob += quality_match_probs.at(int32_t(base_qualities.at(i)));
                }

            } else if (edit.from_length() == edit.to_length() && !edit.sequence().empty()) {

                for (int32_t i = cur_pos; i < cur_pos + edit.from_length(); ++i) {

                    align_path_prob += quality_mismatch_probs.at(int32_t(base_qualities.at(i)));
                }
            
            } else if (edit.from_length() == 0 && edit.to_length() > 0 && !edit.sequence().empty()) {

                align_path_prob += edit.to_length() * indel_prob;

            } else if (edit.from_length() > 0 && edit.to_length() == 0) {

                align_path_prob += edit.from_length() * indel_prob;
            }
        }
    } 

    return align_path_prob;
}

bool operator==(const ReadPathProbs & lhs, const ReadPathProbs & rhs) { 

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

bool operator!=(const ReadPathProbs & lhs, const ReadPathProbs & rhs) { 

    return !(lhs == rhs);
}

bool operator<(const ReadPathProbs & lhs, const ReadPathProbs & rhs) { 

    if (lhs.noise_prob != rhs.noise_prob) {

        return (lhs.noise_prob < rhs.noise_prob);    
    } 

    assert(lhs.read_path_probs.size() == rhs.read_path_probs.size());

    for (size_t i = 0; i < lhs.read_path_probs.size(); ++i) {

        if (lhs.read_path_probs.at(i) != rhs.read_path_probs.at(i)) {

            return (lhs.read_path_probs.at(i) < rhs.read_path_probs.at(i));    
        }         
    }   

    return false;
}

ostream & operator<<(ostream & os, const ReadPathProbs & probs) {

    os << probs.noise_prob;

    for (auto & prob: probs.read_path_probs) {

        os << " " << prob;
    }

    return os;
}

