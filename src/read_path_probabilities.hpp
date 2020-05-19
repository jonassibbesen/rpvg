
#ifndef RPVG_SRC_READPATHPROBABILITIES_HPP
#define RPVG_SRC_READPATHPROBABILITIES_HPP

#include <vector>

#include "vg/io/basic_stream.hpp"
#include "alignment_path.hpp"
#include "paths_index.hpp"
#include "fragment_length_dist.hpp"
#include "path_cluster_estimates.hpp"
#include "utils.hpp"

using namespace std;


class ReadPathProbabilities {

    public: 
    	
    	ReadPathProbabilities(const uint32_t read_count_in, const uint32_t num_paths, const double score_log_base_in, const FragmentLengthDist & fragment_length_dist_in);

        uint32_t readCount() const;
        double noiseProbability() const;
        const vector<double> & probabilities() const;

        void addReadCount(const uint32_t read_count_in);
        void calcReadPathProbabilities(const vector<AlignmentPath> & align_paths, const unordered_map<uint32_t, uint32_t> & clustered_path_index, const vector<PathInfo> & cluster_paths, const bool is_single_end);

        bool mergeIdenticalReadPathProbabilities(const ReadPathProbabilities & probs_2, const double prob_precision);
        vector<pair<double, vector<uint32_t> > > collapsedProbabilities(const double precision) const;

    private:

        uint32_t read_count;
        double noise_prob;
        vector<double> read_path_probs;
        
        const double score_log_base;
        const FragmentLengthDist & fragment_length_dist;
};

bool operator==(const ReadPathProbabilities & lhs, const ReadPathProbabilities & rhs);
bool operator!=(const ReadPathProbabilities & lhs, const ReadPathProbabilities & rhs);
bool operator<(const ReadPathProbabilities & lhs, const ReadPathProbabilities & rhs);

ostream & operator<<(ostream & os, const ReadPathProbabilities & read_path_probs);


#endif
