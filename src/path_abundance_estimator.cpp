
#include <limits>

#include "path_abundance_estimator.hpp"
#include "discrete_sampler.hpp"


PathAbundanceEstimator::PathAbundanceEstimator(const uint32_t max_em_its_in, const double min_abundance_in) : max_em_its(max_em_its_in), min_abundance(min_abundance_in) {}

Abundances PathAbundanceEstimator::inferPathClusterAbundances(const vector<pair<ReadPathProbabilities, uint32_t> > & cluster_probs, const uint32_t num_paths) {

    if (!cluster_probs.empty()) {

        Abundances abundances(num_paths + 1);

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
        removeNoiseAndRenormalizeAbundances(&abundances);

        return abundances;

    } else {

        return Abundances(num_paths, true);
    }
}

void PathAbundanceEstimator::expectationMaximizationEstimator(Abundances * abundances, const Eigen::ColMatrixXd & read_path_probs, const Eigen::RowVectorXui & read_counts) const {

    abundances->read_count = read_counts.sum();
    assert(abundances->read_count > 0);

    Eigen::RowVectorXd prev_read_counts = abundances->expression * abundances->read_count;

    for (size_t i = 0; i < max_em_its; ++i) {

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

MinimumPathAbundanceEstimator::MinimumPathAbundanceEstimator(const uint32_t max_em_its, const double min_abundance) : PathAbundanceEstimator(max_em_its, min_abundance) {}

Abundances MinimumPathAbundanceEstimator::inferPathClusterAbundances(const vector<pair<ReadPathProbabilities, uint32_t> > & cluster_probs, const uint32_t num_paths) {

    if (!cluster_probs.empty()) {

        Eigen::ColVectorXd noise_probs(cluster_probs.size());
        Eigen::RowVectorXui read_counts(cluster_probs.size());

        Eigen::ColMatrixXb read_path_cover(cluster_probs.size(), num_paths);
        Eigen::RowVectorXd path_weights = Eigen::RowVectorXd::Zero(num_paths);

        for (size_t i = 0; i < read_path_cover.rows(); ++i) {

            noise_probs(i) = cluster_probs.at(i).first.noiseProbability();

            if (doubleCompare(noise_probs(i), 1)) {

                read_counts(i) = 0;

            } else {
                
                read_counts(i) = cluster_probs.at(i).second;
            }

            assert(cluster_probs.at(i).first.probabilities().size() == read_path_cover.cols());

            for (size_t j = 0; j < num_paths; ++j) {

                path_weights(j) += log(cluster_probs.at(i).first.probabilities().at(j) + noise_probs(i)) * read_counts(i);

                if (doubleCompare(cluster_probs.at(i).first.probabilities().at(j), 0)) {

                    read_path_cover(i, j) = false;

                } else {
                    
                    read_path_cover(i, j) = true;
                }
            }
        }

        path_weights *= -1;

        vector<uint32_t> min_path_cover = weightedMinimumPathCover(read_path_cover, read_counts, path_weights);

        if (min_path_cover.empty()) {

            return Abundances(num_paths, true);
        }

        Eigen::ColMatrixXd min_path_read_path_probs(cluster_probs.size(), min_path_cover.size());

        for (size_t i = 0; i < min_path_read_path_probs.rows(); ++i) {

            read_counts(i) = cluster_probs.at(i).second;

            for (size_t j = 0; j < min_path_cover.size(); ++j) {

                min_path_read_path_probs(i, j) = cluster_probs.at(i).first.probabilities().at(min_path_cover.at(j));
            }
        }
        
        min_path_read_path_probs = min_path_read_path_probs.array().colwise() / min_path_read_path_probs.rowwise().sum().array();
        min_path_read_path_probs = min_path_read_path_probs.array().colwise() * (1 - noise_probs.array());
        min_path_read_path_probs = min_path_read_path_probs.array().isNaN().select(0, min_path_read_path_probs);

        min_path_read_path_probs.conservativeResize(min_path_read_path_probs.rows(), min_path_read_path_probs.cols() + 1);
        min_path_read_path_probs.col(min_path_read_path_probs.cols() - 1) = noise_probs;

        assert(min_path_read_path_probs.cols() > 1);
        Abundances min_path_abundances(min_path_read_path_probs.cols());
        
        expectationMaximizationEstimator(&min_path_abundances, min_path_read_path_probs, read_counts);

        Abundances abundances(num_paths + 1);

        for (size_t j = 0; j < min_path_cover.size(); j++) {

            abundances.confidence(min_path_cover.at(j)) = min_path_abundances.confidence(j);
            abundances.expression(min_path_cover.at(j)) = min_path_abundances.expression(j);
        }

        assert(min_path_abundances.confidence.cols() == min_path_cover.size() + 1);

        abundances.confidence(min_path_cover.size()) = min_path_abundances.confidence(min_path_cover.size());
        abundances.expression(min_path_cover.size()) = min_path_abundances.expression(min_path_cover.size());  
                  
        removeNoiseAndRenormalizeAbundances(&abundances);
        return abundances;

    } else {

        return Abundances(num_paths, true);
    }
}

vector<uint32_t> MinimumPathAbundanceEstimator::weightedMinimumPathCover(const Eigen::ColMatrixXb & read_path_cover, const Eigen::RowVectorXui & read_counts, const Eigen::RowVectorXd & path_weights) {

    assert(read_path_cover.rows() == read_counts.cols());
    assert(read_path_cover.cols() == path_weights.cols());

    if (read_path_cover.cols() == 1) {

        return vector<uint32_t>({0});
    }

    auto uncovered_read_counts = read_counts;

    vector<uint32_t> min_path_cover;
    min_path_cover.reserve(read_path_cover.cols());

    while (uncovered_read_counts.maxCoeff() > 0) {

        Eigen::RowVectorXd weighted_read_path_cover = (uncovered_read_counts.cast<double>() * read_path_cover.cast<double>()).array() / path_weights.array();
        assert(weighted_read_path_cover.size() == read_path_cover.cols());

        double max_weighted_read_path_cover = weighted_read_path_cover(0);
        uint32_t max_weighted_read_path_cover_idx = 0;

        for (size_t i = 1; i < weighted_read_path_cover.size(); ++i) {

            if (weighted_read_path_cover(i) > max_weighted_read_path_cover) {

                max_weighted_read_path_cover = weighted_read_path_cover(i);
                max_weighted_read_path_cover_idx = i;
            }
        }

        assert(max_weighted_read_path_cover > 0);
        min_path_cover.emplace_back(max_weighted_read_path_cover_idx);

        uncovered_read_counts = (uncovered_read_counts.array() * (!read_path_cover.col(max_weighted_read_path_cover_idx).transpose().array()).cast<uint32_t>()).matrix();
    }

    assert(min_path_cover.size() <= read_path_cover.cols());
    return min_path_cover;
}

