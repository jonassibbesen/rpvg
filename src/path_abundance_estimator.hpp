
#ifndef RPVG_SRC_PATHABUNDANCEESTIMATOR_HPP
#define RPVG_SRC_PATHABUNDANCEESTIMATOR_HPP

#include <vector>
#include <random>

#include <Eigen/Dense>

#include "path_estimator.hpp"
#include "path_cluster_estimates.hpp"
#include "read_path_probabilities.hpp"
#include "utils.hpp"

using namespace std;


class PathAbundanceEstimator : public PathEstimator {

    public:

        PathAbundanceEstimator(const uint32_t max_em_its_in, const double max_rel_em_conv_in, const uint32_t num_gibbs_samples_in, const uint32_t gibbs_thin_its_in, const double prob_precision);
        virtual ~PathAbundanceEstimator() {};

        void estimate(PathClusterEstimates * path_cluster_estimates, const vector<ReadPathProbabilities> & cluster_probs, mt19937 * mt_rng);

    protected: 

        const uint32_t max_em_its;
        const double max_rel_em_conv;

        const uint32_t num_gibbs_samples;
        const uint32_t gibbs_thin_its;

        void EMAbundanceEstimator(PathClusterEstimates * path_cluster_estimates, const Utils::ColMatrixXd & read_path_probs, const Utils::RowVectorXd & read_counts) const;
        void gibbsReadCountSampler(PathClusterEstimates * path_cluster_estimates, const Utils::ColMatrixXd & read_path_probs, const Utils::RowVectorXd & read_counts, const double gamma, mt19937 * mt_rng, const uint32_t num_samples) const;
};

class MinimumPathAbundanceEstimator : public PathAbundanceEstimator {

    public:

        MinimumPathAbundanceEstimator(const uint32_t max_em_its, const double max_rel_em_conv, const uint32_t num_gibbs_samples, const uint32_t gibbs_thin_its, const double prob_precision);
        ~MinimumPathAbundanceEstimator() {};

        void estimate(PathClusterEstimates * path_cluster_estimates, const vector<ReadPathProbabilities> & cluster_probs, mt19937 * mt_rng);

        vector<uint32_t> weightedMinimumPathCover(const Utils::ColMatrixXb & read_path_cover, const Utils::RowVectorXd & read_counts, const Utils::RowVectorXd & path_weights) const;
};

class NestedPathAbundanceEstimator : public PathAbundanceEstimator {

    public:

        NestedPathAbundanceEstimator(const uint32_t group_size_in, const double min_hap_prob_in, const bool infer_collapsed_in, const bool use_group_post_gibbs_in, const uint32_t max_em_its, const double max_rel_em_conv, const uint32_t num_gibbs_samples, const uint32_t gibbs_thin_its, const double prob_precision);
        ~NestedPathAbundanceEstimator() {};

        void estimate(PathClusterEstimates * path_cluster_estimates, const vector<ReadPathProbabilities> & cluster_probs, mt19937 * mt_rng);

    private:

        const uint32_t group_size;
        const double min_hap_prob;

        const bool infer_collapsed;
        const bool use_group_post_gibbs;

        void inferAbundancesIndependentGroups(PathClusterEstimates * path_cluster_estimates, const vector<ReadPathProbabilities> & cluster_probs, mt19937 * mt_rng) const;        
        void inferAbundancesCollapsedGroups(PathClusterEstimates * path_cluster_estimates, const vector<ReadPathProbabilities> & cluster_probs, mt19937 * mt_rng);        

        vector<vector<uint32_t> > findPathGroups(const vector<PathInfo> & paths) const;
        pair<vector<vector<uint32_t> >, vector<uint32_t> > findPathSourceGroups(const vector<PathInfo> & paths) const;

        void sampleGroupPathIndices(vector<vector<uint32_t> > * path_subset_samples, const PathClusterEstimates & group_path_cluster_estimates, const vector<uint32_t> & group, mt19937 * mt_rng) const;
        void selectPathSubsetIndices(spp::sparse_hash_map<vector<uint32_t>, double> * path_subset_samples, const PathClusterEstimates & group_path_cluster_estimates, const vector<vector<uint32_t> > & path_groups, mt19937 * mt_rng) const;

        void inferPathSubsetAbundance(PathClusterEstimates * path_cluster_estimates, const vector<ReadPathProbabilities> & cluster_probs, mt19937 * mt_rng, const spp::sparse_hash_map<vector<uint32_t>, double> & path_subset_samples) const;
};

 
#endif
