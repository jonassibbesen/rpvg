
#include <limits>

#include "path_abundance_estimator.hpp"


PathAbundanceEstimator::PathAbundanceEstimator(const double min_abundance_in, const uint32_t rng_seed) : min_abundance(min_abundance_in) {

    mt_rng = mt19937(rng_seed);
}

void PathAbundanceEstimator::nullifyAbundances(Abundances * abundances) const {

    for (size_t i = 0; i < abundances->expression.cols(); ++i) {

        abundances->confidence(i) = 0;
        abundances->expression(i) = 0;            
    }

    abundances->read_count = 0;    
}

void PathAbundanceEstimator::removeNoiseAndRenormalizeAbundances(Abundances * abundances) const {

    const double noise_read_count = abundances->expression(0, abundances->expression.cols() - 1) * abundances->read_count;
    assert(abundances->read_count >= noise_read_count);

    abundances->confidence.conservativeResize(1, abundances->confidence.cols() - 1);
    abundances->expression.conservativeResize(1, abundances->expression.cols() - 1);

    if (abundances->expression.sum() > 0) {

        abundances->expression = abundances->expression / abundances->expression.sum();
    } 

    abundances->read_count -= noise_read_count;
}

SimplePathAbundanceEstimator::SimplePathAbundanceEstimator(const uint32_t max_em_iterations_in, const double min_abundance, const uint32_t rng_seed) : max_em_iterations(max_em_iterations_in), PathAbundanceEstimator(min_abundance, rng_seed) {}

Abundances SimplePathAbundanceEstimator::inferPathClusterAbundances(const vector<pair<ReadPathProbabilities, uint32_t> > & cluster_probs, const uint32_t num_paths) const {

    Abundances abundances(num_paths + 1);

    if (!cluster_probs.empty()) {

        Eigen::ColMatrixXd read_path_probs(cluster_probs.size(), num_paths + 1);
        Eigen::RowVectorXui read_counts(cluster_probs.size());

        for (size_t i = 0; i < read_path_probs.rows(); ++i) {

            read_counts(i) = cluster_probs.at(i).second;
            assert(cluster_probs.at(i).first.probabilities().size() == read_path_probs.cols() - 1);

            for (size_t j = 0; j < num_paths; ++j) {

                read_path_probs(i, j) = cluster_probs.at(i).first.probabilities().at(j);
            }

            read_path_probs(i, num_paths) = cluster_probs.at(i).first.noiseProbability();
        }

        expectationMaximizationEstimator(&abundances, read_path_probs, read_counts);
    
    } else {

        nullifyAbundances(&abundances);
    }

    removeNoiseAndRenormalizeAbundances(&abundances);
    return abundances;
}

void SimplePathAbundanceEstimator::expectationMaximizationEstimator(Abundances * abundances, const Eigen::ColMatrixXd & read_path_probs, const Eigen::RowVectorXui & read_counts) const {

    abundances->read_count = read_counts.sum();
    assert(abundances->read_count > 0);

    Eigen::RowVectorXd prev_read_counts = abundances->expression * abundances->read_count;

    for (size_t i = 0; i < max_em_iterations; ++i) {

        Eigen::ColMatrixXd posteriors = read_path_probs.array().rowwise() * abundances->expression.array();
        posteriors = posteriors.array().colwise() / posteriors.rowwise().sum().array();

        abundances->expression = read_counts.cast<double>() * posteriors;

        if ((abundances->expression.array() - prev_read_counts.array()).abs().maxCoeff() < min_abundance * abundances->read_count) {

            break;
        } 

        prev_read_counts = abundances->expression;
        abundances->expression /= abundances->read_count;   
    }

    abundances->expression = abundances->expression / abundances->expression.sum();

    for (size_t i = 0; i < abundances->expression.cols(); ++i) {

        if (abundances->expression(i) < min_abundance) {

            abundances->confidence(i) = 0;
            abundances->expression(i) = 0;            
        }
    }

    if (abundances->expression.sum() > 0) {

        abundances->expression = abundances->expression / abundances->expression.sum();
    } 
}

MinimumPathAbundanceEstimator::MinimumPathAbundanceEstimator(const uint32_t num_path_iterations_in, const uint32_t max_em_iterations, const double min_abundance, const uint32_t rng_seed) : num_path_iterations(num_path_iterations_in), SimplePathAbundanceEstimator(max_em_iterations, min_abundance, rng_seed) {}

Abundances MinimumPathAbundanceEstimator::inferPathClusterAbundances(const vector<pair<ReadPathProbabilities, uint32_t> > & cluster_probs, const uint32_t num_paths) {

    Abundances abundances(num_paths + 1);

    if (!cluster_probs.empty()) {

        Eigen::ColMatrixXd read_path_log_probs(cluster_probs.size(), num_paths);
        Eigen::RowVectorXui read_counts(cluster_probs.size());

        for (size_t i = 0; i < read_path_log_probs.rows(); ++i) {

            read_counts(i) = cluster_probs.at(i).second;
            assert(cluster_probs.at(i).first.probabilities().size() == read_path_log_probs.cols());

            for (size_t j = 0; j < num_paths; ++j) {

                read_path_log_probs(i, j) = log(cluster_probs.at(i).first.probabilities().at(j));
            }
        }

        nullifyAbundances(&abundances);

        abundances.read_count = read_counts.sum();
        assert(abundances.read_count > 0);

        for (size_t i = 0; i < num_path_iterations; ++i) {

            Eigen::ColMatrixXd min_path_read_path_probs;

    
    // label_probs = label_probs.array().colwise() / label_probs.rowwise().sum().array();
    // label_probs = label_probs.array().colwise() * (1 - inference_unit.noiseProbability()s.array());
    // label_probs = label_probs.array().isNaN().select(0, label_probs);

    // label_probs.conservativeResize(label_probs.rows(), label_probs.cols() + 1);
    // label_probs.col(label_probs.cols() - 1) = inference_unit.noiseProbability()s;

            vector<uint32_t> min_column_cover = sampleMinimumPathCover(read_path_log_probs, read_counts);

            Abundances min_path_abundances(min_column_cover.size());
            expectationMaximizationEstimator(&min_path_abundances, min_path_read_path_probs, read_counts);

            assert(min_column_cover.size() == min_path_read_path_probs.cols());

            for (size_t j = 0; j < min_column_cover.size(); j++) {

                if (min_path_abundances.confidence(j) > 0) {

                    assert(doubleCompare(min_path_abundances.confidence(j), 1));
                    abundances.confidence(min_column_cover.at(j)) += min_path_abundances.confidence(j);
                    abundances.expression(min_column_cover.at(j)) += min_path_abundances.expression(j);
                }
            }

            assert(min_path_abundances.read_count == abundances.read_count);
        }

        abundances.expression = abundances.expression.array() / abundances.confidence.array();
        abundances.expression = abundances.expression.array().isNaN().select(0, abundances.expression);
        abundances.confidence = abundances.confidence / num_path_iterations;

    } else {

        nullifyAbundances(&abundances);
    }

    removeNoiseAndRenormalizeAbundances(&abundances);
    return abundances;
}

vector<uint32_t> MinimumPathAbundanceEstimator::sampleMinimumPathCover(const Eigen::ColMatrixXd & read_path_log_probs, const Eigen::RowVectorXui & read_counts) {

    auto uncovered_read_counts = read_counts;

    vector<uint32_t> min_path_cover;
    min_path_cover.reserve(read_path_log_probs.cols());

    while (uncovered_read_counts.maxCoeff() > 0) {

        Eigen::RowVectorXd read_path_log_prob_cover = uncovered_read_counts.cast<double>() * read_path_log_probs;
        assert(read_path_log_prob_cover.size() == read_path_log_probs.cols());

        double read_path_log_prob_cover_sum = numeric_limits<double>::lowest();

        for (auto & log_probs: read_path_log_prob_cover) {

            if (doubleCompare(log_probs, 0)) {

                log_probs = numeric_limits<double>::lowest();
            }

            read_path_log_prob_cover_sum = add_log(read_path_log_prob_cover_sum, log_probs);
        }

        vector<double> read_path_prob_cover;
        read_path_prob_cover.reserve(read_path_log_prob_cover.size());

        for (auto & log_probs: read_path_log_prob_cover) {

            read_path_prob_cover.emplace_back(exp(log_probs - read_path_log_prob_cover_sum));
        }

        discrete_distribution<uint32_t> discrete_dist(read_path_prob_cover.begin(), read_path_prob_cover.end());

        const uint32_t sampled_path_idx = discrete_dist(mt_rng); 
        min_path_cover.emplace_back(sampled_path_idx);

        uncovered_read_counts = (uncovered_read_counts.array() * (!read_path_log_probs.col(sampled_path_idx).transpose().cast<bool>().array()).cast<uint32_t>()).matrix();
    }

    return min_path_cover;
}