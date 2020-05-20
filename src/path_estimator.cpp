
#include "path_estimator.hpp"


bool probabilityCountRowsSorter(const pair<Eigen::RowVectorXd, uint32_t> & lhs, const pair<Eigen::RowVectorXd, uint32_t> & rhs) { 

    assert(lhs.first.cols() == rhs.first.cols());

    for (size_t i = 0; i < lhs.first.cols(); ++i) {

        if (!doubleCompare(lhs.first(i), rhs.first(i))) {

            return (lhs.first(i) < rhs.first(i));    
        }         
    }   

    if (lhs.second != rhs.second) {

        return (lhs.second < rhs.second);
    }

    return false;
}

PathEstimator::PathEstimator(const double prob_precision_in) : prob_precision(prob_precision_in) {}

void PathEstimator::constructProbabilityMatrix(Eigen::ColMatrixXd * read_path_probs, Eigen::ColVectorXd * noise_probs, Eigen::RowVectorXui * read_counts, const vector<ReadPathProbabilities> & cluster_probs, const bool add_noise) {

    assert(!cluster_probs.empty());

    *read_path_probs = Eigen::ColMatrixXd(cluster_probs.size(), cluster_probs.front().probabilities().size() + static_cast<uint32_t>(add_noise));
    *noise_probs = Eigen::ColVectorXd(cluster_probs.size());
    *read_counts = Eigen::RowVectorXui(cluster_probs.size());

    for (size_t i = 0; i < read_path_probs->rows(); ++i) {

        assert(cluster_probs.at(i).probabilities().size() + static_cast<uint32_t>(add_noise) == read_path_probs->cols());

        for (size_t j = 0; j < cluster_probs.at(i).probabilities().size(); ++j) {

            (*read_path_probs)(i, j) = cluster_probs.at(i).probabilities().at(j);
        }

        if (add_noise) {

            (*read_path_probs)(i, cluster_probs.at(i).probabilities().size()) = cluster_probs.at(i).noiseProbability();
        }

        (*noise_probs)(i) = cluster_probs.at(i).noiseProbability();
        (*read_counts)(i) = cluster_probs.at(i).readCount();
    } 
}

void PathEstimator::addNoiseAndNormalizeProbabilityMatrix(Eigen::ColMatrixXd * read_path_probs, const Eigen::ColVectorXd & noise_probs) {

    assert(read_path_probs->rows() == noise_probs.rows());

    *read_path_probs = read_path_probs->array().colwise() / read_path_probs->rowwise().sum().array();
    *read_path_probs = read_path_probs->array().colwise() * (1 - noise_probs.array());
    *read_path_probs = read_path_probs->array().isNaN().select(0, *read_path_probs);

    read_path_probs->conservativeResize(read_path_probs->rows(), read_path_probs->cols() + 1);
    read_path_probs->col(read_path_probs->cols() - 1) = noise_probs;
}

void PathEstimator::sortProbabilityMatrix(Eigen::ColMatrixXd * read_path_probs, Eigen::RowVectorXui * read_counts) {

    assert(read_path_probs->rows() > 0);
    assert(read_path_probs->rows() == read_counts->cols());

    vector<pair<Eigen::RowVectorXd, uint32_t> > read_path_prob_rows;
    read_path_prob_rows.reserve(read_path_probs->rows());

    for (size_t i = 0; i < read_path_probs->rows(); ++i) {

        read_path_prob_rows.emplace_back(read_path_probs->row(i), (*read_counts)(i));
    }

    sort(read_path_prob_rows.begin(), read_path_prob_rows.end(), probabilityCountRowsSorter);

    for (size_t i = 0; i < read_path_probs->rows(); ++i) {
    
        read_path_probs->row(i) = read_path_prob_rows.at(i).first;
        (*read_counts)(i) = read_path_prob_rows.at(i).second;
    }    
}

void PathEstimator::collapseProbabilityMatrix(Eigen::ColMatrixXd * read_path_probs, Eigen::RowVectorXui * read_counts) {

    assert(read_path_probs->rows() > 0);
    assert(read_path_probs->rows() == read_counts->cols());

    uint32_t prev_unique_probs_row = 0;

    for (size_t i = 1; i < read_path_probs->rows(); ++i) {

        bool is_identical = true;

        for (size_t j = 0; j < read_path_probs->cols(); ++j) {

            if (abs((*read_path_probs)(prev_unique_probs_row, j) - (*read_path_probs)(i, j)) >= prob_precision) {

                is_identical = false;
                break;
            }
        }

        if (is_identical) {

            (*read_counts)(prev_unique_probs_row) += (*read_counts)(i);

        } else {

            if (prev_unique_probs_row + 1 < i) {

                read_path_probs->row(prev_unique_probs_row + 1) = read_path_probs->row(i);
                (*read_counts)(prev_unique_probs_row + 1) = (*read_counts)(i);
            }

            prev_unique_probs_row++;
        }
    }

    read_path_probs->conservativeResize(prev_unique_probs_row + 1, read_path_probs->cols());
    read_counts->conservativeResize(read_counts->rows(), prev_unique_probs_row + 1);
}
