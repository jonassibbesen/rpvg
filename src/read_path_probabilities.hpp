
#ifndef RPVG_SRC_READPATHPROBABILITIES_HPP
#define RPVG_SRC_READPATHPROBABILITIES_HPP

#include <vector>

#include "sparsepp/spp.h"

#include "vg/io/basic_stream.hpp"
#include "alignment_path.hpp"
#include "paths_index.hpp"
#include "fragment_length_dist.hpp"
#include "path_cluster_estimates.hpp"
#include "utils.hpp"

using namespace std;


class ReadPathProbabilities {

    public: 

        ReadPathProbabilities();
    	ReadPathProbabilities(const uint32_t read_count_in, const double prob_precision_in);

        uint32_t readCount() const;
        double noiseProb() const;
        const vector<pair<double, vector<uint32_t> > > & pathProbs() const;

        static vector<double> calcAlignPathLogProbs(const vector<AlignmentPath> & align_paths, const FragmentLengthDist & fragment_length_dist, const bool is_single_end);

        void addReadCount(const uint32_t read_count_in);
        void addPathProbs(const vector<AlignmentPath> & align_paths, const vector<vector<gbwt::size_type> > & align_paths_ids, const spp::sparse_hash_map<uint32_t, uint32_t> & clustered_path_index, const vector<PathInfo> & cluster_paths, const FragmentLengthDist & fragment_length_dist, const bool is_single_end, const double min_noise_prob, const bool collapse_groups = false, const spp::sparse_hash_map<string, uint32_t> & group_name_index = spp::sparse_hash_map<string, uint32_t>());

        bool quickMergeIdentical(const ReadPathProbabilities & probs_2);

    private:

        uint32_t read_count;
        double noise_prob;
        vector<pair<double, vector<uint32_t> > > path_probs;
        
        double prob_precision;
};

bool operator==(const ReadPathProbabilities & lhs, const ReadPathProbabilities & rhs);
bool operator!=(const ReadPathProbabilities & lhs, const ReadPathProbabilities & rhs);
bool operator<(const ReadPathProbabilities & lhs, const ReadPathProbabilities & rhs);

ostream & operator<<(ostream & os, const ReadPathProbabilities & read_path_probs);


#endif


