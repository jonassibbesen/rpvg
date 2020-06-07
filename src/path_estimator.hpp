
#ifndef RPVG_SRC_PATHESTIMATOR_HPP
#define RPVG_SRC_PATHESTIMATOR_HPP

#include <vector>

#include <Eigen/Dense>

#include "path_cluster_estimates.hpp"
#include "read_path_probabilities.hpp"
#include "utils.hpp"

using namespace std;


class PathEstimator {

    public:

        PathEstimator(const double prob_precision_in);
        virtual ~PathEstimator() {};

        virtual void estimate(PathClusterEstimates * path_cluster_estimates, const vector<ReadPathProbabilities> & cluster_probs) = 0;

    protected:
       
        const double prob_precision;
 
        void constructProbabilityMatrix(Eigen::ColMatrixXd * read_path_probs, Eigen::ColVectorXd * noise_probs, Eigen::RowVectorXui * read_counts, const vector<ReadPathProbabilities> & cluster_probs, const bool add_noise, const double max_noise_prob);
        void addNoiseAndNormalizeProbabilityMatrix(Eigen::ColMatrixXd * read_path_probs, const Eigen::ColVectorXd & noise_probs);

        void collapseProbabilityMatrixReads(Eigen::ColMatrixXd * read_path_probs, Eigen::RowVectorXui * read_counts);
        void collapseProbabilityMatrixPaths(Eigen::ColMatrixXd * read_path_probs);

        void calculatePathGroupPosteriors(PathClusterEstimates * path_cluster_estimates, const Eigen::ColMatrixXd & read_path_probs, const uint32_t noise_col_idx, const Eigen::RowVectorXui & read_counts, const uint32_t group_size);
        void samplePathGroupPosteriorsGibbs(PathClusterEstimates * path_cluster_estimates, const Eigen::ColMatrixXd & read_path_probs, const uint32_t noise_col_idx, const Eigen::RowVectorXui & read_counts, const uint32_t group_size, const uint32_t num_gibbs_its, mt19937 * mt_rng);

    private:

        void rowSortProbabilityMatrix(Eigen::ColMatrixXd * read_path_probs, Eigen::RowVectorXui * read_counts);
        void colSortProbabilityMatrix(Eigen::ColMatrixXd * read_path_probs);
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
