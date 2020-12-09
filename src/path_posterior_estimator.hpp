
#ifndef RPVG_SRC_PATHLIKELIHOODESTIMATOR_HPP
#define RPVG_SRC_PATHLIKELIHOODESTIMATOR_HPP

#include <vector>
#include <random>

#include <Eigen/Dense>

#include "path_estimator.hpp"
#include "path_cluster_estimates.hpp"
#include "read_path_probabilities.hpp"
#include "utils.hpp"

using namespace std;


class PathPosteriorEstimator : public PathEstimator {

    public:

        PathPosteriorEstimator(const double prob_precision);
        virtual ~PathPosteriorEstimator() {};

        void estimate(PathClusterEstimates * path_cluster_estimates, const vector<ReadPathProbabilities> & cluster_probs, mt19937 * mt_rng);

};

class PathGroupPosteriorEstimator : public PathPosteriorEstimator {

    public:

        PathGroupPosteriorEstimator(const uint32_t group_size_in, const bool use_group_post_gibbs_in, const double prob_precision);
        ~PathGroupPosteriorEstimator() {};

        void estimate(PathClusterEstimates * path_cluster_estimates, const vector<ReadPathProbabilities> & cluster_probs, mt19937 * mt_rng);

    private: 

        const uint32_t group_size;
        const bool use_group_post_gibbs;
};

 
#endif
