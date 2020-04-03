
#include "path_abundance_estimator.hpp"


PathAbundanceEstimator::PathAbundanceEstimator(const double min_abundance_in) : min_abundance(min_abundance_in) {}

void PathAbundanceEstimator::removeNoiseAndRenormalize(Abundances * abundances) const {

    const double noise_read_count = abundances->expression(0, abundances->expression.cols() - 1) * abundances->read_count;
    assert(abundances->read_count >= noise_read_count);

    abundances->confidence.conservativeResize(1, abundances->confidence.cols() - 1);
    abundances->expression.conservativeResize(1, abundances->expression.cols() - 1);

    if (abundances->expression.sum() > 0) {

        abundances->expression = abundances->expression / abundances->expression.sum();
    } 

    abundances->read_count -= noise_read_count;
}

EMPathAbundanceEstimator::EMPathAbundanceEstimator(const double min_abundance_in, const uint32_t max_em_iteration_in) : PathAbundanceEstimator(min_abundance_in), max_em_iteration(max_em_iteration_in) {}

Abundances EMPathAbundanceEstimator::inferPathClusterAbundances(const vector<pair<ReadPathProbabilities, uint32_t> > & cluster_probs, const uint32_t num_paths) const {

    Abundances abundances(num_paths + 1);

    if (!cluster_probs.empty()) {

        Eigen::RowVectorXi read_counts(cluster_probs.size());
        Eigen::RowMatrixXd read_path_probs(cluster_probs.size(), num_paths + 1);

        for (size_t i = 0; i < read_path_probs.rows(); ++i) {

            read_counts(i) = cluster_probs.at(i).second;
            assert(cluster_probs.at(i).first.read_path_probs.size() == read_path_probs.cols() - 1);

            for (size_t j = 0; j < num_paths; ++j) {

                read_path_probs(i, j) = cluster_probs.at(i).first.read_path_probs.at(j);
            }

            read_path_probs(i, num_paths) = cluster_probs.at(i).first.noise_prob;
        }

        read_path_probs = read_path_probs.array().colwise() / read_path_probs.rowwise().sum().array();

        abundances.read_count = read_counts.sum();
        assert(abundances.read_count > 0);

        Eigen::RowVectorXd prev_read_counts = abundances.expression * abundances.read_count;

        for (size_t i = 0; i < max_em_iteration; ++i) {

            Eigen::RowMatrixXd posteriors = read_path_probs.array().rowwise() * abundances.expression.array();
            posteriors = posteriors.array().colwise() / posteriors.rowwise().sum().array();

            abundances.expression = read_counts.cast<double>() * posteriors;

            if ((abundances.expression.array() - prev_read_counts.array()).abs().maxCoeff() < min_abundance * abundances.read_count) {

                break;
            } 

            prev_read_counts = abundances.expression;
            abundances.expression /= abundances.read_count;   
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
    
    } else {

        for (size_t i = 0; i < abundances.expression.cols(); ++i) {

            abundances.confidence(i) = 0;
            abundances.expression(i) = 0;            
        }

        abundances.read_count = 0;
    }

    removeNoiseAndRenormalize(&abundances);
    return abundances;
}

