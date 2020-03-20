
#ifndef RPVG_SRC_READPATHPROBS_HPP
#define RPVG_SRC_READPATHPROBS_HPP

#include <vector>

#include "vg/io/basic_stream.hpp"
#include "alignment_path.hpp"
#include "paths_index.hpp"
#include "fragment_length_dist.hpp"
#include "utils.hpp"


using namespace std;


class ReadPathProbs {

    public: 
    	
    	ReadPathProbs();
    	ReadPathProbs(const int32_t num_paths, const double score_log_base_in);
        
        double noise_prob;
        vector<double> read_path_probs;

        void calcReadPathProbs(const vector<AlignmentPath> & align_paths, const unordered_map<int32_t, int32_t> & clustered_path_index, const FragmentLengthDist & fragment_length_dist, const bool is_single_end);
        void addPositionalProbs(const vector<double> & path_lengths);

    private:

        double score_log_base;

    	double calcReadMappingProbs(const vg::Alignment & alignment, const vector<double> & quality_match_probs, const vector<double> & quality_mismatch_probs, const double indel_prob) const;
};

bool operator==(const ReadPathProbs & lhs, const ReadPathProbs & rhs);
bool operator!=(const ReadPathProbs & lhs, const ReadPathProbs & rhs);
bool operator<(const ReadPathProbs & lhs, const ReadPathProbs & rhs);

ostream & operator<<(ostream & os, const ReadPathProbs & probs);


#endif
