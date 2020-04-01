
#include "abundance_estimator.hpp"


AbundanceEstimator::AbundanceEstimator(const double min_abundance_in) : min_abundance(min_abundance_in) {}

EMAbundanceEstimator::EMAbundanceEstimator(const double min_abundance_in, const uint32_t max_em_iteration_in, const double stop_em_count_diff_in) : AbundanceEstimator(min_abundance_in), max_em_iteration(max_em_iteration_in), stop_em_count_diff(stop_em_count_diff_in) {}

Abundances EMAbundanceEstimator::inferClusterAbundance(const vector<pair<ReadPathProbabilities, uint32_t> > & cluster_probs, const uint32_t num_paths) {

    Abundances abundances(num_paths + 1);

    if (!cluster_probs.empty()) {

        Eigen::RowVectorXi read_counts(cluster_probs.size());
        Eigen::RowMatrixXd read_path_probs(cluster_probs.size(), num_paths + 1);

        for (size_t i = 0; i < read_path_probs.rows(); ++i) {

            read_counts(i) = cluster_probs.at(i).second;
            assert(cluster_probs.at(i).first.read_path_probs.size() == read_path_probs.cols() - 1);

            for (size_t j = 0; j < read_path_probs.cols() - 1; ++j) {

                read_path_probs(i, j) = cluster_probs.at(i).first.read_path_probs.at(j);
            }

            read_path_probs(i, read_path_probs.cols() - 1) = cluster_probs.at(i).first.noise_prob;
        }

        read_path_probs = read_path_probs.array().colwise() / read_path_probs.rowwise().sum().array();

        const uint32_t total_read_count = read_counts.sum();
        Eigen::RowVectorXd prev_read_counts = abundances.expression * total_read_count;

        for (size_t i = 0; i < max_em_iteration; ++i) {

            Eigen::RowMatrixXd posteriors = read_path_probs.array().rowwise() * abundances.expression.array();
            posteriors = posteriors.array().colwise() / posteriors.rowwise().sum().array();

            abundances.expression = read_counts.cast<double>() * posteriors;

            if ((abundances.expression.array() - prev_read_counts.array()).abs().maxCoeff() < stop_em_count_diff) {

                break;
            } 

            prev_read_counts = abundances.expression;
            abundances.expression /= total_read_count;   
        }

        abundances.expression = abundances.expression / abundances.expression.sum();

        for (size_t i = 0; i < abundances.expression.cols(); ++i) {

            if (abundances.expression(i) < min_abundance) {

                abundances.confidence(i) = 0;
                abundances.expression(i) = 0;            
            }
        }

        if (abundances.expression.sum() > 0) {

            abundances.expression = abundances.expression / abundances.expression.sum();
        } 

        return abundances;
    
    } else {

        for (size_t i = 0; i < abundances.expression.cols(); ++i) {

            if (abundances.expression(i) < min_abundance) {

                abundances.confidence(i) = 0;
                abundances.expression(i) = 0;            
            }
        }
    }

    return abundances;
}

