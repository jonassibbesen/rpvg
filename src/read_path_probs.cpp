
#include "read_path_probs.hpp"
#include "utils.hpp"


ReadPathProbs::ReadPathProbs() {

    noise_prob = 1;    
}

ReadPathProbs::ReadPathProbs(const int32_t num_paths) {

    noise_prob = 1;
    read_path_probs = vector<double>(num_paths, 0); 
}

void ReadPathProbs::calcReadPathProbs(const vector<AlignmentPath> & align_paths, const unordered_map<uint32_t, uint32_t> & clustered_path_index, const uint32_t frag_length_mean, const uint32_t frag_length_sd) {

    assert(!align_paths.empty());
    assert(clustered_path_index.size() == read_path_probs.size());

    if (align_paths.front().mapqs.first > 0 && align_paths.front().mapqs.second > 0) {

        noise_prob = mapqsToProb(align_paths.front().mapqs);
        assert(noise_prob < 1);

        vector<double> align_paths_probs;
        align_paths_probs.reserve(align_paths.size());

        if (align_paths.size() == 1) {

            align_paths_probs.emplace_back(1 - noise_prob);
        
        } else {

            double sum_scores = 0;

            for (auto & align_path: align_paths) {

                sum_scores += align_path.scores.first + align_path.scores.second;
            }

            double sum_probs = 0;

            for (auto & align_path: align_paths) {

                align_paths_probs.emplace_back((align_path.scores.first + align_path.scores.second) / sum_scores * normal_pdf<double>(align_path.seq_length, frag_length_mean, frag_length_sd));
                sum_probs += align_paths_probs.back();
            }

            for (auto & probs: align_paths_probs) {

                probs /= sum_probs;
                probs *= (1 - noise_prob);
            }
        }

        for (size_t i = 0; i < align_paths.size(); ++i) {

            for (auto & path: align_paths.at(i).path_ids) {

                read_path_probs.at(clustered_path_index.at(path)) = align_paths_probs.at(i) / align_paths.at(i).path_ids.size();
            }
        }
    }
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

ostream& operator<<(ostream& os, const ReadPathProbs & probs) {

    os << probs.noise_prob;

    for (auto & prob: probs.read_path_probs) {

        os << " " << prob;
    }

    os << endl;

    return os;
}


