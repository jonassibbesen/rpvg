
#ifndef RPVG_SRC_READPATHPROBABILITIES_HPP
#define RPVG_SRC_READPATHPROBABILITIES_HPP

#include <vector>

#include "vg/io/basic_stream.hpp"
#include "alignment_path.hpp"
#include "paths_index.hpp"
#include "fragment_length_dist.hpp"
#include "path.hpp"
#include "utils.hpp"

using namespace std;


class ReadPathProbabilities {

    public: 
    	
    	ReadPathProbabilities();
    	ReadPathProbabilities(const uint32_t num_paths, const double score_log_base_in, const FragmentLengthDist & fragment_length_dist_in);
        
        const vector<double> & probabilities() const;
        double noiseProbability() const;

        void calcReadPathProbabilities(const vector<AlignmentPath> & align_paths, const unordered_map<uint32_t, uint32_t> & clustered_path_index, const vector<Path> & cluster_paths, const bool is_single_end);

        vector<pair<double, vector<uint32_t> > > collapsedProbabilities(const double precision) const;

    private:

        double noise_prob;
        vector<double> read_path_probs;

        const double score_log_base;
        const FragmentLengthDist fragment_length_dist;
};

bool operator==(const ReadPathProbabilities & lhs, const ReadPathProbabilities & rhs);
bool operator!=(const ReadPathProbabilities & lhs, const ReadPathProbabilities & rhs);
bool operator<(const ReadPathProbabilities & lhs, const ReadPathProbabilities & rhs);

ostream & operator<<(ostream & os, const ReadPathProbabilities & read_path_probs);


#endif
