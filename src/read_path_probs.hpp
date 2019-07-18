
#ifndef VGPROB_SRC_READPATHPROBS_HPP
#define VGPROB_SRC_READPATHPROBS_HPP

#include <vector>

#include "vg/io/basic_stream.hpp"
#include "alignment_path.hpp"
#include "fragment_length_dist.hpp"
#include "utils.hpp"


using namespace std;


class ReadPathProbs {

    public: 
    	
    	ReadPathProbs();
    	ReadPathProbs(const int32_t num_paths);
        
        void calcReadPathProbs(const vector<AlignmentPath> & align_paths, const unordered_map<uint32_t, uint32_t> & clustered_path_index, const FragmentLengthDist & fragment_length_dist);

        double score_log_base;
        double noise_prob;
        vector<double> read_path_probs;

    private:

    	double calcReadMappingProbs(const vg::Alignment & alignment, const vector<double> & quality_match_probs, const vector<double> & quality_mismatch_probs, const double indel_prob) const;
};

inline bool operator==(const ReadPathProbs & lhs, const ReadPathProbs & rhs) { 

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

inline bool operator!=(const ReadPathProbs & lhs, const ReadPathProbs & rhs) { 

    return !(lhs == rhs);
}

inline bool operator<(const ReadPathProbs & lhs, const ReadPathProbs & rhs) { 

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

inline ostream & operator<<(ostream & os, const ReadPathProbs & probs) {

    os << probs.noise_prob;

    for (auto & prob: probs.read_path_probs) {

        os << " " << prob;
    }

    return os;
}


#endif
