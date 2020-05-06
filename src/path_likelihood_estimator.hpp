
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


class PathLikelihoodEstimator : public PathEstimator {

    public:

        PathLikelihoodEstimator(const bool use_log_in, const double prob_precision);
        ~PathLikelihoodEstimator() {};

        void estimate(PathClusterEstimates * path_cluster_estimates, const vector<ReadPathProbabilities> & cluster_probs);

    protected: 

        const bool use_log; 

};

class PathGroupLikelihoodEstimator : public PathLikelihoodEstimator {

    public:

        PathGroupLikelihoodEstimator(const uint32_t ploidy_in, const bool use_log, const double prob_precision);
        ~PathGroupLikelihoodEstimator() {};

        void estimate(PathClusterEstimates * path_cluster_estimates, const vector<ReadPathProbabilities> & cluster_probs);

    private: 

        const uint32_t ploidy;
};

 
#endif
