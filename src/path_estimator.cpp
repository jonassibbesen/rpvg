
#include "path_estimator.hpp"

static const double max_log_likelihood_diff = log(pow(10, -16));

static const uint32_t min_gibbs_chains = 10;
static const double gibbs_chain_scaling = 0.01;

static const uint32_t min_burn_it = 50;
static const double burn_it_scaling = 0.025;

static const uint32_t min_gibbs_it = 100; 
static const double gibbs_it_scaling = 0.05; 

bool probabilityCountRowSorter(const pair<Utils::RowVectorXd, uint32_t> & lhs, const pair<Utils::RowVectorXd, uint32_t> & rhs) { 

    assert(lhs.first.cols() == rhs.first.cols());

    for (size_t i = 0; i < lhs.first.cols(); ++i) {

        if (!Utils::doubleCompare(lhs.first(i), rhs.first(i))) {

            return (lhs.first(i) < rhs.first(i));    
        }         
    }   

    if (lhs.second != rhs.second) {

        return (lhs.second < rhs.second);
    }

    return false;
}

bool probabilityCountColSorter(const pair<Utils::ColVectorXd, uint32_t> & lhs, const pair<Utils::ColVectorXd, uint32_t> & rhs) { 

    assert(lhs.first.rows() == rhs.first.rows());

    for (size_t i = 0; i < lhs.first.rows(); ++i) {

        if (!Utils::doubleCompare(lhs.first(i), rhs.first(i))) {

            return (lhs.first(i) < rhs.first(i));    
        }         
    }   

    if (lhs.second != rhs.second) {

        return (lhs.second < rhs.second);
    }

    return false;
}

PathEstimator::PathEstimator(const double prob_precision_in) : prob_precision(prob_precision_in) {}

void PathEstimator::constructProbabilityMatrix(Utils::ColMatrixXd * read_path_probs, Utils::ColVectorXd * noise_probs, Utils::RowVectorXui * read_counts, const vector<ReadPathProbabilities> & cluster_probs, const uint32_t num_paths) const {

    assert(!cluster_probs.empty());

    *read_path_probs = Utils::ColMatrixXd::Zero(cluster_probs.size(), num_paths);
    *noise_probs = Utils::ColVectorXd(cluster_probs.size());
    *read_counts = Utils::RowVectorXui(cluster_probs.size());

    for (size_t i = 0; i < cluster_probs.size(); ++i) {

        for (auto & path_probs: cluster_probs.at(i).pathProbs()) {

            for (auto & path: path_probs.second) {

                assert(path < num_paths);
                (*read_path_probs)(i, path) = path_probs.first;
            }
        }

        (*noise_probs)(i, 0) = cluster_probs.at(i).noiseProb();
        (*read_counts)(0, i) = cluster_probs.at(i).readCount();
    }
}

void PathEstimator::constructPartialProbabilityMatrix(Utils::ColMatrixXd * read_path_probs, Utils::ColVectorXd * noise_probs, Utils::RowVectorXui * read_counts, const vector<ReadPathProbabilities> & cluster_probs, const vector<uint32_t> & path_ids, const uint32_t num_paths) const {

    assert(!cluster_probs.empty());
    assert(!path_ids.empty());

    vector<int32_t> path_id_idx(num_paths, -1);

    for (size_t i = 0; i < path_ids.size(); ++i) {

        path_id_idx.at(path_ids.at(i)) = i;
    }

    *read_path_probs = Utils::ColMatrixXd::Zero(cluster_probs.size(), path_ids.size());
    *noise_probs = Utils::ColVectorXd(cluster_probs.size());
    *read_counts = Utils::RowVectorXui(cluster_probs.size());

    for (size_t i = 0; i < cluster_probs.size(); ++i) {

        for (auto & path_probs: cluster_probs.at(i).pathProbs()) {

            for (auto & path: path_probs.second) {
    
                assert(path < num_paths);

                if (path_id_idx.at(path) >= 0) {

                    (*read_path_probs)(i, path_id_idx.at(path)) = path_probs.first;
                }
            }
        }

        (*noise_probs)(i, 0) = cluster_probs.at(i).noiseProb();
        (*read_counts)(0, i) = cluster_probs.at(i).readCount();
    }
}

void PathEstimator::constructGroupedProbabilityMatrix(Utils::ColMatrixXd * read_path_probs, Utils::ColVectorXd * noise_probs, Utils::RowVectorXui * read_counts, const vector<ReadPathProbabilities> & cluster_probs, const vector<vector<uint32_t> > & path_groups, const uint32_t num_paths) const {

    assert(!cluster_probs.empty());
    assert(!path_groups.empty());

    vector<vector<uint32_t> > path_id_group_idx(num_paths, vector<uint32_t>());

    for (size_t i = 0; i < path_groups.size(); ++i) {

        assert(!path_groups.at(i).empty());

        for (auto & path: path_groups.at(i)) {    

            path_id_group_idx.at(path).emplace_back(i);
        }
    }

    *read_path_probs = Utils::ColMatrixXd::Zero(cluster_probs.size(), path_groups.size());
    *noise_probs = Utils::ColVectorXd(cluster_probs.size());
    *read_counts = Utils::RowVectorXui(cluster_probs.size());

    for (size_t i = 0; i < cluster_probs.size(); ++i) {

        for (auto & path_probs: cluster_probs.at(i).pathProbs()) {

            for (auto & path: path_probs.second) {

                assert(path < num_paths);

                for (auto & group_id: path_id_group_idx.at(path)) {

                    (*read_path_probs)(i, group_id) += path_probs.first;
                }
            }
        }

        (*noise_probs)(i, 0) = cluster_probs.at(i).noiseProb();
        (*read_counts)(0, i) = cluster_probs.at(i).readCount();
    }
}

void PathEstimator::addNoiseAndNormalizeProbabilityMatrix(Utils::ColMatrixXd * read_path_probs, const Utils::ColVectorXd & noise_probs) const {

    assert(read_path_probs->rows() == noise_probs.rows());

    *read_path_probs = read_path_probs->array().colwise() / read_path_probs->rowwise().sum().array();
    *read_path_probs = read_path_probs->array().colwise() * (1 - noise_probs.array());
    *read_path_probs = read_path_probs->array().isNaN().select(0, *read_path_probs);

    read_path_probs->conservativeResize(read_path_probs->rows(), read_path_probs->cols() + 1);
    read_path_probs->col(read_path_probs->cols() - 1) = noise_probs;
}

void PathEstimator::rowSortProbabilityMatrix(Utils::ColMatrixXd * read_path_probs, Utils::RowVectorXui * read_counts) const {

    assert(read_path_probs->rows() > 0);
    assert(read_path_probs->rows() == read_counts->cols());

    vector<pair<Utils::RowVectorXd, uint32_t> > read_path_prob_rows;
    read_path_prob_rows.reserve(read_path_probs->rows());

    for (size_t i = 0; i < read_path_probs->rows(); ++i) {

        read_path_prob_rows.emplace_back(read_path_probs->row(i), (*read_counts)(0, i));
    }

    sort(read_path_prob_rows.begin(), read_path_prob_rows.end(), probabilityCountRowSorter);

    for (size_t i = 0; i < read_path_probs->rows(); ++i) {
    
        read_path_probs->row(i) = read_path_prob_rows.at(i).first;
        (*read_counts)(0, i) = read_path_prob_rows.at(i).second;
    }    
}

void PathEstimator::readCollapseProbabilityMatrix(Utils::ColMatrixXd * read_path_probs, Utils::RowVectorXui * read_counts) const {

    assert(read_path_probs->rows() > 0);
    assert(read_path_probs->rows() == read_counts->cols());

    rowSortProbabilityMatrix(read_path_probs, read_counts);

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

            read_counts->col(prev_unique_probs_row) += read_counts->col(i);

        } else {

            if (prev_unique_probs_row + 1 < i) {

                read_path_probs->row(prev_unique_probs_row + 1) = read_path_probs->row(i);
                read_counts->col(prev_unique_probs_row + 1) = read_counts->col(i);
            }

            prev_unique_probs_row++;
        }
    }

    read_path_probs->conservativeResize(prev_unique_probs_row + 1, read_path_probs->cols());
    read_counts->conservativeResize(read_counts->rows(), prev_unique_probs_row + 1);
}

void PathEstimator::colSortProbabilityMatrix(Utils::ColMatrixXd * read_path_probs) const {

    assert(read_path_probs->cols() > 0);

    vector<pair<Utils::ColVectorXd, uint32_t> > read_path_prob_cols;
    read_path_prob_cols.reserve(read_path_probs->cols());

    for (size_t i = 0; i < read_path_probs->cols(); ++i) {

        read_path_prob_cols.emplace_back(read_path_probs->col(i), i);
    }

    sort(read_path_prob_cols.begin(), read_path_prob_cols.end(), probabilityCountColSorter);

    for (size_t i = 0; i < read_path_probs->cols(); ++i) {
    
        read_path_probs->col(i) = read_path_prob_cols.at(i).first;
    }    
}

void PathEstimator::pathCollapseProbabilityMatrix(Utils::ColMatrixXd * read_path_probs) const {

    assert(read_path_probs->cols() > 0);    
    colSortProbabilityMatrix(read_path_probs);

    uint32_t prev_unique_probs_col = 0;

    for (size_t i = 1; i < read_path_probs->cols(); ++i) {

        bool is_identical = true;

        for (size_t j = 0; j < read_path_probs->rows(); ++j) {

            if (abs((*read_path_probs)(j, prev_unique_probs_col) - (*read_path_probs)(j, i)) >= prob_precision) {

                is_identical = false;
                break;
            }
        }

        if (!is_identical) {

            if (prev_unique_probs_col + 1 < i) {

                read_path_probs->col(prev_unique_probs_col + 1) = read_path_probs->col(i);
            }

            prev_unique_probs_col++;
        }
    }

    read_path_probs->conservativeResize(read_path_probs->rows(), prev_unique_probs_col + 1);
}

vector<double> PathEstimator::calcPathLogFrequences(const vector<uint32_t> & path_counts) const {

    vector<double> path_log_freqs;
    path_log_freqs.reserve(path_counts.size());

    uint32_t count_sum = accumulate(path_counts.begin(), path_counts.end(), 0);
    assert(count_sum > 0);

    for (auto & count: path_counts) {

        assert(count > 0);
        path_log_freqs.emplace_back(log(count / static_cast<double>(count_sum)));
    }

    return path_log_freqs;
}

void PathEstimator::calculatePathGroupPosteriorsFull(PathClusterEstimates * path_cluster_estimates, const Utils::ColMatrixXd & read_path_probs, const Utils::ColVectorXd & noise_probs, const Utils::RowVectorXui & read_counts, const vector<uint32_t> & path_counts, const uint32_t group_size) const {

    assert(read_path_probs.rows() > 0);
    assert(read_path_probs.rows() == noise_probs.rows());
    assert(read_path_probs.rows() == read_counts.cols());
    assert(read_path_probs.cols() == path_counts.size());
    assert(group_size > 0);

    auto path_log_freqs = calcPathLogFrequences(path_counts);
    assert(path_log_freqs.size() == path_counts.size());

    path_cluster_estimates->initEstimates(read_path_probs.cols(), group_size, true);

    assert(path_cluster_estimates->posteriors.cols() > 0);
    assert(path_cluster_estimates->posteriors.cols() == path_cluster_estimates->path_group_sets.size());

    double sum_log_posterior = numeric_limits<double>::lowest();

    for (uint32_t i = 0; i < path_cluster_estimates->path_group_sets.size(); ++i) {

        assert(path_cluster_estimates->path_group_sets.at(i).size() == group_size);

        Utils::ColVectorXd group_read_probs = noise_probs;

        for (auto & path_idx: path_cluster_estimates->path_group_sets.at(i)) {

            group_read_probs += (read_path_probs.col(path_idx) / static_cast<double>(group_size));
        }

        path_cluster_estimates->posteriors(0, i) = read_counts.cast<double>() * group_read_probs.array().log().matrix();

        for (auto & path_idx: path_cluster_estimates->path_group_sets.at(i)) {
            
            path_cluster_estimates->posteriors(0, i) += path_log_freqs.at(path_idx);
        }

        path_cluster_estimates->posteriors(0, i) += log(Utils::numPermutations(path_cluster_estimates->path_group_sets.at(i)));

        sum_log_posterior = Utils::add_log(sum_log_posterior, path_cluster_estimates->posteriors(0, i));
    }

    for (size_t i = 0; i < path_cluster_estimates->posteriors.cols(); ++i) {

        path_cluster_estimates->posteriors(0, i) = exp(path_cluster_estimates->posteriors(0, i) - sum_log_posterior);
    }
}

void PathEstimator::calculatePathGroupPosteriorsBounded(PathClusterEstimates * path_cluster_estimates, const Utils::ColMatrixXd & read_path_probs, const Utils::ColVectorXd & noise_probs, const Utils::RowVectorXui & read_counts, const vector<uint32_t> & path_counts, const uint32_t group_size) const {

    assert(read_path_probs.rows() > 0);
    assert(read_path_probs.rows() == noise_probs.rows());
    assert(read_path_probs.rows() == read_counts.cols());
    assert(read_path_probs.cols() == path_counts.size());
    assert(group_size == 2);

    auto path_log_freqs = calcPathLogFrequences(path_counts);
    assert(path_log_freqs.size() == path_counts.size());

    path_cluster_estimates->initEstimates(0, 0, true);

    assert(path_cluster_estimates->posteriors.cols() == 0);
    assert(path_cluster_estimates->posteriors.cols() == path_cluster_estimates->path_group_sets.size());

    PathClusterEstimates marginal_path_cluster_estimates;
    calculatePathGroupPosteriorsFull(&marginal_path_cluster_estimates, read_path_probs, noise_probs, read_counts, path_counts, 1);

    assert(marginal_path_cluster_estimates.posteriors.cols() == read_path_probs.cols());
    assert(marginal_path_cluster_estimates.posteriors.cols() == marginal_path_cluster_estimates.path_group_sets.size());

    vector<pair<double, uint32_t> > marginal_posteriors;
    marginal_posteriors.reserve(marginal_path_cluster_estimates.posteriors.cols());

    for (size_t i = 0; i < marginal_path_cluster_estimates.posteriors.cols(); ++i) {

        assert(marginal_path_cluster_estimates.path_group_sets.at(i).size() == 1);
        marginal_posteriors.emplace_back(marginal_path_cluster_estimates.posteriors(0, i), marginal_path_cluster_estimates.path_group_sets.at(i).front());
    }

    sort(marginal_posteriors.rbegin(), marginal_posteriors.rend());

    const Utils::ColVectorXd max_read_probs = (read_path_probs.rowwise().maxCoeff() / static_cast<double>(group_size));

    vector<double> log_likelihoods;

    double max_log_likelihood = numeric_limits<double>::lowest(); 
    double sum_log_posterior = numeric_limits<double>::lowest();

    for (uint32_t i = 0; i < marginal_posteriors.size(); ++i) {

        const uint32_t first_path_idx = marginal_posteriors.at(i).second;

        Utils::ColVectorXd group_read_probs_base = noise_probs;
        group_read_probs_base += (read_path_probs.col(first_path_idx) / static_cast<double>(group_size));

        double max_log_likelihood_best = read_counts.cast<double>() * (group_read_probs_base + max_read_probs).array().log().matrix();
        max_log_likelihood_best += path_log_freqs.at(first_path_idx) + log(2);

        if (max_log_likelihood_best - max_log_likelihood < max_log_likelihood_diff) {

            continue;
        }

        for (uint32_t j = i; j < marginal_posteriors.size(); ++j) {

            const uint32_t second_path_idx = marginal_posteriors.at(j).second;

            log_likelihoods.emplace_back(read_counts.cast<double>() * (group_read_probs_base + (read_path_probs.col(second_path_idx) / static_cast<double>(group_size))).array().log().matrix());
            log_likelihoods.back() += path_log_freqs.at(first_path_idx) + path_log_freqs.at(second_path_idx) + log(Utils::numPermutations(vector<uint32_t>({first_path_idx, second_path_idx})));

            if (log_likelihoods.back() - max_log_likelihood < max_log_likelihood_diff) {

                log_likelihoods.pop_back();
                continue;
            }

            max_log_likelihood = max(max_log_likelihood, log_likelihoods.back());
            sum_log_posterior = Utils::add_log(sum_log_posterior, log_likelihoods.back());

            path_cluster_estimates->path_group_sets.emplace_back(vector<uint32_t>({first_path_idx, second_path_idx}));
        }
    }

    assert(log_likelihoods.size() == path_cluster_estimates->path_group_sets.size());
    path_cluster_estimates->posteriors = Utils::RowVectorXd(1, log_likelihoods.size());

    for (size_t i = 0; i < path_cluster_estimates->posteriors.cols(); ++i) {

        path_cluster_estimates->posteriors(0, i) = exp(log_likelihoods.at(i) - sum_log_posterior);
    }
}

void PathEstimator::estimatePathGroupPosteriorsGibbs(PathClusterEstimates * path_cluster_estimates, const Utils::ColMatrixXd & read_path_probs, const Utils::ColVectorXd & noise_probs, const Utils::RowVectorXui & read_counts, const vector<uint32_t> & path_counts, const uint32_t group_size, mt19937 * mt_rng) const {

    assert(read_path_probs.rows() > 0);
    assert(read_path_probs.rows() == noise_probs.rows());
    assert(read_path_probs.rows() == read_counts.cols());
    assert(read_path_probs.cols() == path_counts.size());
    assert(group_size > 0);

    auto path_log_freqs = calcPathLogFrequences(path_counts);
    assert(path_log_freqs.size() == path_counts.size());

    path_cluster_estimates->initEstimates(0, 0, true);

    assert(path_cluster_estimates->posteriors.cols() == 0);
    assert(path_cluster_estimates->posteriors.cols() == path_cluster_estimates->path_group_sets.size());

    uniform_int_distribution<uint32_t> init_path_sampler(0, path_log_freqs.size() - 1);

    vector<uint32_t> cur_sampled_group_paths;
    cur_sampled_group_paths.reserve(group_size);

    spp::sparse_hash_map<vector<uint32_t>, discrete_distribution<uint32_t> > group_path_sampler_cache;
    spp::sparse_hash_map<vector<uint32_t>, uint32_t> path_group_sets_indices;

    vector<uint32_t> path_group_sample_counts;

    const uint32_t num_gibbs_chains = min_gibbs_chains + round(gibbs_chain_scaling * group_size * path_log_freqs.size());
    const uint32_t num_burn_its = min_burn_it + round(burn_it_scaling * group_size * path_log_freqs.size());
    const uint32_t num_gibbs_its = min_gibbs_it + round(gibbs_it_scaling * group_size * path_log_freqs.size());

    for (uint32_t c = 0; c < num_gibbs_chains; ++c) {

        cur_sampled_group_paths.clear();

        for (uint32_t i = 0; i < group_size; ++i) {

            cur_sampled_group_paths.emplace_back(init_path_sampler(*mt_rng));
        }

        for (uint32_t i = 0; i < num_burn_its + num_gibbs_its; ++i) {

            for (uint32_t j = 0; j < group_size; ++j) {

                vector<uint32_t> new_cur_sampled_group_paths = cur_sampled_group_paths;

                new_cur_sampled_group_paths.at(j) = read_path_probs.cols();
                sort(new_cur_sampled_group_paths.begin(), new_cur_sampled_group_paths.end());

                auto group_path_sampler_cache_it = group_path_sampler_cache.emplace(new_cur_sampled_group_paths, discrete_distribution<uint32_t>());

                if (group_path_sampler_cache_it.second) {

                    Utils::ColVectorXd group_read_probs = noise_probs;

                    for (uint32_t k = 0; k < group_size; ++k) {

                        if (j != k) {

                            group_read_probs += (read_path_probs.col(cur_sampled_group_paths.at(k)) / static_cast<double>(group_size));
                        }
                    }

                    vector<double> group_probs;
                    group_probs.reserve(read_path_probs.cols());

                    double sum_log_group_probs = numeric_limits<double>::lowest();

                    for (uint32_t k = 0; k < read_path_probs.cols(); ++k) {

                        group_probs.emplace_back(read_counts.cast<double>() * (group_read_probs + (read_path_probs.col(k) / static_cast<double>(group_size))).array().log().matrix());
                        group_probs.back() += path_log_freqs.at(k);

                        sum_log_group_probs = Utils::add_log(sum_log_group_probs, group_probs.back());
                    }

                    for (auto & prob: group_probs) {

                        prob = exp(prob - sum_log_group_probs);
                    }

                    group_path_sampler_cache_it.first->second = discrete_distribution<uint32_t>(group_probs.begin(), group_probs.end());
                }

                cur_sampled_group_paths.at(j) = group_path_sampler_cache_it.first->second(*mt_rng);
            }

            if (i >= num_burn_its) {

                vector<uint32_t> cur_sampled_group_paths_sort = cur_sampled_group_paths;
                sort(cur_sampled_group_paths_sort.begin(), cur_sampled_group_paths_sort.end());

                auto path_group_sets_indices_it = path_group_sets_indices.emplace(cur_sampled_group_paths_sort, path_cluster_estimates->path_group_sets.size());

                if (path_group_sets_indices_it.second) {

                    path_cluster_estimates->path_group_sets.emplace_back(cur_sampled_group_paths_sort);
                    path_group_sample_counts.emplace_back(1);

                } else {

                    path_group_sample_counts.at(path_group_sets_indices_it.first->second)++;
                }
            }
        }
    }

    path_cluster_estimates->posteriors = Utils::RowVectorXd(1, path_group_sample_counts.size());

    for (size_t i = 0; i < path_cluster_estimates->posteriors.cols(); ++i) {

        path_cluster_estimates->posteriors(0, i) = path_group_sample_counts.at(i) / static_cast<double>(num_gibbs_chains * num_gibbs_its);
    }
}

