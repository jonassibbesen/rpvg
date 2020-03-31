
#ifndef RPVG_SRC_READPATHPROBABILITIES_HPP
#define RPVG_SRC_READPATHPROBABILITIES_HPP

#include <vector>

#include "vg/io/basic_stream.hpp"
#include "alignment_path.hpp"
#include "paths_index.hpp"
#include "fragment_length_dist.hpp"
#include "utils.hpp"

using namespace std;


class ReadPathProbabilities {

    public: 
    	
    	ReadPathProbabilities();
    	ReadPathProbabilities(const uint32_t num_paths, const double score_log_base_in);
        
        double noise_prob;

        vector<double> read_path_probs;
        vector<double> positive_prob_indices;

        void calcReadPathProbabilities(const vector<AlignmentPath> & align_paths, const unordered_map<uint32_t, uint32_t> & clustered_path_index, const FragmentLengthDist & fragment_length_dist, const bool is_single_end);
        void addPositionalProbabilities(const vector<double> & path_lengths);

    private:

        double score_log_base;

    	double calcReadMappingProbabilities(const vg::Alignment & alignment, const vector<double> & quality_match_probs, const vector<double> & quality_mismatch_probs, const double indel_prob) const;
};

bool operator==(const ReadPathProbabilities & lhs, const ReadPathProbabilities & rhs);
bool operator!=(const ReadPathProbabilities & lhs, const ReadPathProbabilities & rhs);
bool operator<(const ReadPathProbabilities & lhs, const ReadPathProbabilities & rhs);

ostream & operator<<(ostream & os, const ReadPathProbabilities & probs);


#endif
