
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

        virtual void estimate(PathClusterEstimates * path_cluster_estimates, const vector<ReadPathProbabilities> & cluster_probs, mt19937 * mt_rng) = 0;

    protected:
       
        const double prob_precision;

        void constructProbabilityMatrix(Utils::ColMatrixXd * read_path_probs, Utils::ColVectorXd * noise_probs, Utils::RowVectorXd * read_counts, const vector<ReadPathProbabilities> & cluster_probs, const uint32_t num_paths) const; 
        void constructPartialProbabilityMatrix(Utils::ColMatrixXd * read_path_probs, Utils::ColVectorXd * noise_probs, Utils::RowVectorXd * read_counts, const vector<ReadPathProbabilities> & cluster_probs, const vector<uint32_t> & path_ids, const uint32_t num_paths) const;
        void constructGroupedProbabilityMatrix(Utils::ColMatrixXd * read_path_probs, Utils::ColVectorXd * noise_probs, Utils::RowVectorXd * read_counts, const vector<ReadPathProbabilities> & cluster_probs, const vector<vector<uint32_t> > & path_groups, const uint32_t num_paths) const;

        void addNoiseAndNormalizeProbabilityMatrix(Utils::ColMatrixXd * read_path_probs, const Utils::ColVectorXd & noise_probs) const;
        void detractNoiseAndNormalizeProbabilityMatrix(Utils::ColMatrixXd * read_path_probs, Utils::ColVectorXd * noise_probs, Utils::RowVectorXd * read_counts) const;

        void readCollapseProbabilityMatrix(Utils::ColMatrixXd * read_path_probs, Utils::RowVectorXd * read_counts) const;
        void pathCollapseProbabilityMatrix(Utils::ColMatrixXd * read_path_probs) const;

        void calculatePathGroupPosteriorsFull(PathClusterEstimates * path_cluster_estimates, const Utils::ColMatrixXd & read_path_probs, const Utils::ColVectorXd & noise_probs, const Utils::RowVectorXd & read_counts, const vector<uint32_t> & path_counts, const uint32_t group_size) const;
        void calculatePathGroupPosteriorsBounded(PathClusterEstimates * path_cluster_estimates, const Utils::ColMatrixXd & read_path_probs, const Utils::ColVectorXd & noise_probs, const Utils::RowVectorXd & read_counts, const vector<uint32_t> & path_counts, const uint32_t group_size, const double min_rel_likelihood) const;
        void estimatePathGroupPosteriorsGibbs(PathClusterEstimates * path_cluster_estimates, const Utils::ColMatrixXd & read_path_probs, const Utils::ColVectorXd & noise_probs, const Utils::RowVectorXd & read_counts, const vector<uint32_t> & path_counts, const uint32_t group_size, mt19937 * mt_rng) const;

    private:

        void rowSortProbabilityMatrix(Utils::ColMatrixXd * read_path_probs, Utils::RowVectorXd * read_counts) const;
        void colSortProbabilityMatrix(Utils::ColMatrixXd * read_path_probs) const;

        vector<double> calcPathLogFrequences(const vector<uint32_t> & path_counts) const;
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
