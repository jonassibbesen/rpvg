
#ifndef RPVG_SRC_PATHLIKELIHOODESTIMATOR_HPP
#define RPVG_SRC_PATHLIKELIHOODESTIMATOR_HPP

#include <vector>
#include <random>

#include <Eigen/Dense>

#include "path_likelihoods.hpp"
#include "read_path_probabilities.hpp"
#include "utils.hpp"

using namespace std;


class PathLikelihoodEstimator {

    public:

        PathLikelihoodEstimator(const double prob_precision_in);
        ~PathLikelihoodEstimator() {};

        const double prob_precision;

        PathLikelihoods inferPathClusterLogLikelihoods(const vector<ReadPathProbabilities> & cluster_probs, const vector<PathInfo> & cluster_paths);

        static void constructProbabilityMatrix(Eigen::ColMatrixXd * read_path_probs, Eigen::ColVectorXd * noise_probs, Eigen::RowVectorXui * read_counts, const vector<ReadPathProbabilities> & cluster_probs);
        static void addNoiseToProbabilityMatrix(Eigen::ColMatrixXd * read_path_probs, const Eigen::ColVectorXd & noise_probs);

        static void sortProbabilityMatrix(Eigen::ColMatrixXd * read_path_probs, Eigen::RowVectorXui * read_counts);
        static void collapseProbabilityMatrix(Eigen::ColMatrixXd * read_path_probs, Eigen::RowVectorXui * read_counts, const double collapse_prob_precision);
};

class PathComboLikelihoodEstimator : PathLikelihoodEstimator {

    public:

        PathComboLikelihoodEstimator(const uint32_t ploidy_in, const double prob_precision);
        ~PathComboLikelihoodEstimator() {};

        const uint32_t ploidy;

        PathComboLikelihoods inferPathClusterComboLogLikelihoods(const vector<ReadPathProbabilities> & cluster_probs, const vector<PathInfo> & cluster_paths);
};

namespace std {

    template<> 
    struct hash<vector<uint32_t> >
    {
        size_t operator()(vector<uint32_t> const & vec) const
        {
            size_t seed = 0;

            for (auto & val: vec) {

                spp::hash_combine(seed, val);
            }

            return seed;
        }
    };
}

 
#endif
