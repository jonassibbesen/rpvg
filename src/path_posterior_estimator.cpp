
#include "path_posterior_estimator.hpp"


PathPosteriorEstimator::PathPosteriorEstimator(const double prob_precision) : PathEstimator(prob_precision) {}

void PathPosteriorEstimator::estimate(PathClusterEstimates * path_cluster_estimates, const vector<ReadPathProbabilities> & cluster_probs) {

    if (!cluster_probs.empty()) {

        Eigen::ColMatrixXd read_path_probs;
        Eigen::ColVectorXd noise_probs;
        Eigen::RowVectorXui read_counts;

        constructProbabilityMatrix(&read_path_probs, &noise_probs, &read_counts, cluster_probs, true, 2);

        noise_probs = read_path_probs.col(read_path_probs.cols() - 1);
        read_path_probs.conservativeResize(read_path_probs.rows(), read_path_probs.cols() - 1);

        calculatePathGroupPosteriors(path_cluster_estimates, read_path_probs, noise_probs, read_counts, 1);

        assert(path_cluster_estimates->posteriors.cols() == read_path_probs.cols());
        assert(path_cluster_estimates->posteriors.cols() == path_cluster_estimates->paths.size());
        assert(path_cluster_estimates->posteriors.cols() == path_cluster_estimates->path_groups.size());

    } else {

        path_cluster_estimates->initEstimates(path_cluster_estimates->paths.size(), 1, true);
    }
}

PathGroupPosteriorEstimator::PathGroupPosteriorEstimator(const uint32_t ploidy_in, const bool use_exact_in, const uint32_t rng_seed, const double prob_precision) : ploidy(ploidy_in), use_exact(use_exact_in), PathPosteriorEstimator(prob_precision) {

    mt_rng = mt19937(rng_seed);
}

void PathGroupPosteriorEstimator::estimate(PathClusterEstimates * path_cluster_estimates, const vector<ReadPathProbabilities> & cluster_probs) {

    if (!cluster_probs.empty()) {

        Eigen::ColMatrixXd read_path_probs;
        Eigen::ColVectorXd noise_probs;
        Eigen::RowVectorXui read_counts;

        constructProbabilityMatrix(&read_path_probs, &noise_probs, &read_counts, cluster_probs, true, 2);

        noise_probs = read_path_probs.col(read_path_probs.cols() - 1);
        read_path_probs.conservativeResize(read_path_probs.rows(), read_path_probs.cols() - 1);

        if (use_exact) {

            calculatePathGroupPosteriors(path_cluster_estimates, read_path_probs, noise_probs, read_counts, ploidy);
        
        } else {

            estimatePathGroupPosteriorsGibbs(path_cluster_estimates, read_path_probs, noise_probs, read_counts, ploidy, &mt_rng);
        }

        assert(path_cluster_estimates->posteriors.cols() == path_cluster_estimates->path_groups.size());

    } else {

        path_cluster_estimates->initEstimates(path_cluster_estimates->paths.size(), ploidy, true);
    }
}
