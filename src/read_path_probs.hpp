
#ifndef RPVG_SRC_READPATHPROBS_HPP
#define RPVG_SRC_READPATHPROBS_HPP

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
        
        double noise_prob;
        vector<double> read_path_probs;

        void calcReadPathProbs(const vector<AlignmentPath> & align_paths, const unordered_map<int32_t, int32_t> & clustered_path_index, const FragmentLengthDist & fragment_length_dist);

    private:

        double score_log_base;

    	double calcReadMappingProbs(const vg::Alignment & alignment, const vector<double> & quality_match_probs, const vector<double> & quality_mismatch_probs, const double indel_prob) const;
        vector<double> calcRelativeAlignmentScoreLogProbs(const vector<AlignmentPath> & align_paths) const;
};

bool operator==(const ReadPathProbs & lhs, const ReadPathProbs & rhs);
bool operator!=(const ReadPathProbs & lhs, const ReadPathProbs & rhs);
bool operator<(const ReadPathProbs & lhs, const ReadPathProbs & rhs);

ostream & operator<<(ostream & os, const ReadPathProbs & probs);


#endif
