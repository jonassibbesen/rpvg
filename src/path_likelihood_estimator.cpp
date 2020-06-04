
#include "path_likelihood_estimator.hpp"


PathLikelihoodEstimator::PathLikelihoodEstimator(const bool use_log_in, const double prob_precision) : use_log(use_log_in), PathEstimator(prob_precision) {}

void PathLikelihoodEstimator::estimate(PathClusterEstimates * path_cluster_estimates, const vector<ReadPathProbabilities> & cluster_probs) {

    if (!cluster_probs.empty()) {

        Eigen::ColMatrixXd read_path_probs;
        Eigen::ColVectorXd noise_probs;
        Eigen::RowVectorXui read_counts;

        constructProbabilityMatrix(&read_path_probs, &noise_probs, &read_counts, cluster_probs, true);
        rowCollapseProbabilityMatrix(&read_path_probs, &read_counts);

        path_cluster_estimates->likelihoods = Likelihoods(path_cluster_estimates->paths.size(), 1, use_log);
        assert(path_cluster_estimates->likelihoods.likelihoods.cols() + 1 == read_path_probs.cols());

        for (size_t i = 0; i < path_cluster_estimates->likelihoods.likelihoods.cols(); ++i) {

            path_cluster_estimates->likelihoods.likelihoods(0, i) = read_counts.cast<double>() * (read_path_probs.col(i) + read_path_probs.col(read_path_probs.cols() - 1)).array().log().matrix();
        }

    } else {

        path_cluster_estimates->likelihoods = Likelihoods(path_cluster_estimates->paths.size(), 1, use_log);
    }
}

PathGroupLikelihoodEstimator::PathGroupLikelihoodEstimator(const uint32_t ploidy_in, const bool use_log, const double prob_precision) : ploidy(ploidy_in), PathLikelihoodEstimator(use_log, prob_precision) {

    assert(ploidy >= 1 && ploidy <= 2);
}

void PathGroupLikelihoodEstimator::estimate(PathClusterEstimates * path_cluster_estimates, const vector<ReadPathProbabilities> & cluster_probs) {

    if (!cluster_probs.empty()) {

        if (ploidy == 1) {

            PathLikelihoodEstimator::estimate(path_cluster_estimates, cluster_probs);
            return;
        }

        Eigen::ColMatrixXd read_path_probs;
        Eigen::ColVectorXd noise_probs;
        Eigen::RowVectorXui read_counts;

        constructProbabilityMatrix(&read_path_probs, &noise_probs, &read_counts, cluster_probs, true);
        rowCollapseProbabilityMatrix(&read_path_probs, &read_counts);

        path_cluster_estimates->likelihoods = Likelihoods(path_cluster_estimates->paths.size(), ploidy, use_log);
        assert(path_cluster_estimates->likelihoods.likelihoods.cols() == path_cluster_estimates->likelihoods.groups.size());

        for (uint32_t i = 0; i < path_cluster_estimates->likelihoods.groups.size(); ++i) {

            Eigen::ColVectorXd path_group_probs = read_path_probs.col(read_path_probs.cols() - 1);

            for (auto path_idx: path_cluster_estimates->likelihoods.groups.at(i)) {

                path_group_probs += read_path_probs.col(path_idx);
            }

            path_cluster_estimates->likelihoods.likelihoods(0, i) = read_counts.cast<double>() * path_group_probs.array().log().matrix();
        }

    } else {

        path_cluster_estimates->likelihoods = Likelihoods(path_cluster_estimates->paths.size(), ploidy, use_log);
    }
}

