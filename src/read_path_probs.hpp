
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

        vector<double> calcRelativeAlignmentScoreLogProbs(const vector<AlignmentPath> & align_path) const;
    	double calcReadMappingProbs(const vg::Alignment & alignment, const vector<double> & quality_match_probs, const vector<double> & quality_mismatch_probs, const double indel_prob) const;
};

bool operator==(const ReadPathProbs & lhs, const ReadPathProbs & rhs);
bool operator!=(const ReadPathProbs & lhs, const ReadPathProbs & rhs);
bool operator<(const ReadPathProbs & lhs, const ReadPathProbs & rhs);

ostream & operator<<(ostream & os, const ReadPathProbs & probs);


#endif
