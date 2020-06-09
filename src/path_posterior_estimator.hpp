
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

        void estimate(PathClusterEstimates * path_cluster_estimates, const vector<ReadPathProbabilities> & cluster_probs);

};

class PathGroupPosteriorEstimator : public PathPosteriorEstimator {

    public:

        PathGroupPosteriorEstimator(const uint32_t num_gibbs_its_in, const uint32_t ploidy_in, const bool use_exact_in, const uint32_t rng_seed, const double prob_precision);
        ~PathGroupPosteriorEstimator() {};

        void estimate(PathClusterEstimates * path_cluster_estimates, const vector<ReadPathProbabilities> & cluster_probs);

    private: 

        const uint32_t num_gibbs_its;
        const uint32_t ploidy;
        const bool use_exact;

        mt19937 mt_rng;
};

 
#endif
