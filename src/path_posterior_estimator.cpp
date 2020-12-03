
#include "path_posterior_estimator.hpp"


PathPosteriorEstimator::PathPosteriorEstimator(const double prob_precision) : PathEstimator(prob_precision) {}

void PathPosteriorEstimator::estimate(PathClusterEstimates * path_cluster_estimates, const vector<ReadPathProbabilities> & cluster_probs, mt19937 * mt_rng) {

    if (!cluster_probs.empty()) {

        Eigen::ColMatrixXd read_path_probs;
        Eigen::ColVectorXd noise_probs;
        Eigen::RowVectorXui read_counts;

        constructProbabilityMatrix(&read_path_probs, &noise_probs, &read_counts, cluster_probs, path_cluster_estimates->paths.size());

        vector<uint32_t> path_counts;
        path_counts.reserve(path_cluster_estimates->paths.size());

        for (auto & path: path_cluster_estimates->paths) {

            path_counts.emplace_back(path.source_count);
        }

        calculatePathGroupPosteriorsFull(path_cluster_estimates, read_path_probs, noise_probs, read_counts, path_counts, 1);

        assert(path_cluster_estimates->posteriors.cols() == read_path_probs.cols());
        assert(path_cluster_estimates->posteriors.cols() == path_cluster_estimates->paths.size());
        assert(path_cluster_estimates->posteriors.cols() == path_cluster_estimates->path_group_sets.size());

    } else {

        path_cluster_estimates->initEstimates(path_cluster_estimates->paths.size(), 1, true);
    }
}

PathGroupPosteriorEstimator::PathGroupPosteriorEstimator(const uint32_t ploidy_in, const bool use_hap_gibbs_in, const double prob_precision) : ploidy(ploidy_in), use_hap_gibbs(use_hap_gibbs_in), PathPosteriorEstimator(prob_precision) {}

void PathGroupPosteriorEstimator::estimate(PathClusterEstimates * path_cluster_estimates, const vector<ReadPathProbabilities> & cluster_probs, mt19937 * mt_rng) {

    if (!cluster_probs.empty()) {

        Eigen::ColMatrixXd read_path_probs;
        Eigen::ColVectorXd noise_probs;
        Eigen::RowVectorXui read_counts;

        constructProbabilityMatrix(&read_path_probs, &noise_probs, &read_counts, cluster_probs, path_cluster_estimates->paths.size());

        vector<uint32_t> path_counts;
        path_counts.reserve(path_cluster_estimates->paths.size());

        for (auto & path: path_cluster_estimates->paths) {

            path_counts.emplace_back(path.source_count);
        }

        if (use_hap_gibbs) {

            estimatePathGroupPosteriorsGibbs(path_cluster_estimates, read_path_probs, noise_probs, read_counts, path_counts, ploidy, mt_rng);            

        } else {

            if (ploidy == 2) {

                calculatePathGroupPosteriorsBounded(path_cluster_estimates, read_path_probs, noise_probs, read_counts, path_counts, ploidy);
            
            } else {

                calculatePathGroupPosteriorsFull(path_cluster_estimates, read_path_probs, noise_probs, read_counts, path_counts, ploidy);
            }
        }

        assert(path_cluster_estimates->posteriors.cols() == path_cluster_estimates->path_group_sets.size());

    } else {

        path_cluster_estimates->initEstimates(0, 0, true);
    }
}
