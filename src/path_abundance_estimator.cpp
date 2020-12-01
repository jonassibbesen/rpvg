
#include <limits>
#include <chrono>

#include "path_abundance_estimator.hpp"


const uint32_t min_em_conv_its = 10;
const double min_em_abundances = pow(10, -8);

const double abundance_gibbs_gamma = 1;

PathAbundanceEstimator::PathAbundanceEstimator(const uint32_t max_em_its_in, const double min_em_conv, const uint32_t num_gibbs_samples_in, const uint32_t gibbs_thin_its_in, const double prob_precision) : max_em_its(max_em_its_in), em_conv_min_exp(min_em_conv), em_conv_max_rel_diff(min_em_conv), num_gibbs_samples(num_gibbs_samples_in), gibbs_thin_its(gibbs_thin_its_in), PathEstimator(prob_precision) {}

void PathAbundanceEstimator::estimate(PathClusterEstimates * path_cluster_estimates, const vector<ReadPathProbabilities> & cluster_probs, mt19937 * mt_rng) {

    if (!cluster_probs.empty()) {

        Eigen::ColMatrixXd read_path_probs;
        Eigen::ColVectorXd noise_probs;
        Eigen::RowVectorXui read_counts;

        vector<uint32_t> path_ids(path_cluster_estimates->paths.size());
        iota(path_ids.begin(), path_ids.end(), 0);

        constructProbabilityMatrix(&read_path_probs, &noise_probs, &read_counts, cluster_probs, path_ids);

        read_path_probs.conservativeResize(read_path_probs.rows(), read_path_probs.cols() + 1);
        read_path_probs.col(read_path_probs.cols() - 1) = noise_probs;

        path_cluster_estimates->initEstimates(path_cluster_estimates->paths.size() + 1, 0, false);
        path_cluster_estimates->total_read_count = read_counts.sum();

        EMAbundanceEstimator(path_cluster_estimates, read_path_probs, read_counts);

        if (num_gibbs_samples > 0) {

            vector<CountSamples> * gibbs_read_count_samples = &(path_cluster_estimates->gibbs_read_count_samples);
            gibbs_read_count_samples->emplace_back(CountSamples());

            gibbs_read_count_samples->back().path_ids = vector<uint32_t>(path_cluster_estimates->abundances.cols());
            iota(gibbs_read_count_samples->back().path_ids.begin(), gibbs_read_count_samples->back().path_ids.end(), 0);

            gibbs_read_count_samples->back().samples = vector<vector<double> >(path_cluster_estimates->abundances.cols(), vector<double>());

            gibbsReadCountSampler(path_cluster_estimates, read_path_probs, read_counts, abundance_gibbs_gamma, mt_rng);
        }

        removeNoiseAndRenormalizeAbundances(path_cluster_estimates);

    } else {

        path_cluster_estimates->initEstimates(path_cluster_estimates->paths.size(), 0, true);
    }
}

void PathAbundanceEstimator::EMAbundanceEstimator(PathClusterEstimates * path_cluster_estimates, const Eigen::ColMatrixXd & read_path_probs, const Eigen::RowVectorXui & read_counts) const {

    assert(path_cluster_estimates->total_read_count > 0);

    Eigen::RowVectorXd prev_abundances = path_cluster_estimates->abundances;
    uint32_t em_conv_its = 0;

    for (uint32_t i = 0; i < max_em_its; ++i) {

        Eigen::ColMatrixXd read_posteriors = read_path_probs.array().rowwise() * path_cluster_estimates->abundances.array();
        read_posteriors = read_posteriors.array().colwise() / read_posteriors.rowwise().sum().array();

        path_cluster_estimates->abundances = read_counts.cast<double>() * read_posteriors;
        path_cluster_estimates->abundances /= path_cluster_estimates->total_read_count;

        bool has_converged = true;

        for (size_t i = 0; i < path_cluster_estimates->abundances.cols(); ++i) {

            if (path_cluster_estimates->abundances(0, i) > em_conv_min_exp) {

                auto relative_abundances_diff = fabs(path_cluster_estimates->abundances(0, i) - prev_abundances(0, i)) / path_cluster_estimates->abundances(0, i);

                if (relative_abundances_diff > em_conv_max_rel_diff) {

                    has_converged = false;
                    break;
                }
            }
        }

        if (has_converged) {

            em_conv_its++;

            if (em_conv_its == min_em_conv_its) {

                break;
            }
        
        } else {

            em_conv_its = 0;
        } 

        prev_abundances = path_cluster_estimates->abundances;
    }

    double abundances_sum = 0;

    for (size_t i = 0; i < path_cluster_estimates->abundances.cols(); ++i) {

        if (path_cluster_estimates->abundances(0, i) < min_em_abundances) {

            path_cluster_estimates->abundances(0, i) = 0;                    
        } 

        abundances_sum += path_cluster_estimates->abundances(0, i);
    }

    if (abundances_sum > 0) {

        path_cluster_estimates->abundances = path_cluster_estimates->abundances / abundances_sum;
    }
}

void PathAbundanceEstimator::gibbsReadCountSampler(PathClusterEstimates * path_cluster_estimates, const Eigen::ColMatrixXd & read_path_probs, const Eigen::RowVectorXui & read_counts, const double gamma, mt19937 * mt_rng) const {

    assert(path_cluster_estimates->total_read_count > 0);

    assert(!path_cluster_estimates->gibbs_read_count_samples.empty());
    assert(path_cluster_estimates->gibbs_read_count_samples.back().path_ids.size() == path_cluster_estimates->abundances.cols());
    assert(path_cluster_estimates->gibbs_read_count_samples.back().samples.size() == path_cluster_estimates->abundances.cols());

    assert(doubleCompare(path_cluster_estimates->abundances.sum(), 1));
    Eigen::RowVectorXd gibbs_abundances = path_cluster_estimates->abundances;

    const uint32_t num_gibbs_its = num_gibbs_samples * gibbs_thin_its;

    for (uint32_t gibbs_it = 1; gibbs_it <= num_gibbs_its; ++gibbs_it) {

        Eigen::ColMatrixXd read_posteriors = read_path_probs.array().rowwise() * gibbs_abundances.array();
        read_posteriors = read_posteriors.array().colwise() / read_posteriors.rowwise().sum().array();

        vector<uint32_t> gibbs_path_read_counts(gibbs_abundances.cols(), 0);

        for (size_t i = 0; i < read_path_probs.rows(); ++i) {

            uint32_t row_reads_counts = read_counts(0, i);
            double row_sum_probs = 1;

            for (size_t j = 0; j < read_path_probs.cols(); ++j) {

                auto cur_prob = read_path_probs(i, j);

                if (cur_prob > 0) {

                    assert(row_sum_probs > 0);
                    assert(cur_prob < row_sum_probs || doubleCompare(cur_prob, row_sum_probs));

                    binomial_distribution<uint32_t> path_read_count_sampler(row_reads_counts, max(1.0, cur_prob / row_sum_probs));
                    auto path_read_count = path_read_count_sampler(*mt_rng);

                    gibbs_path_read_counts.at(j) += path_read_count;
                    row_reads_counts -= path_read_count;

                    if (row_reads_counts == 0) {

                        break;
                    }
                }

                row_sum_probs -= cur_prob;
            }

            assert(row_reads_counts == 0);
        }

        double gibbs_abundances_sum = 0;

        for (size_t i = 0; i < gibbs_abundances.cols(); ++i) {

            gamma_distribution<double> gamma_count_dist(gibbs_path_read_counts.at(i) + gamma, 1);

            gibbs_abundances(0, i) = gamma_count_dist(*mt_rng);
            gibbs_abundances_sum += gibbs_abundances(0, i);
        }

        gibbs_abundances = gibbs_abundances / gibbs_abundances_sum;

        if (gibbs_it % gibbs_thin_its == 0) {

            for (size_t i = 0; i < gibbs_abundances.cols(); ++i) {

                path_cluster_estimates->gibbs_read_count_samples.back().samples.at(i).emplace_back(gibbs_abundances(0, i) * path_cluster_estimates->total_read_count);
            }
        }
    }
}

void PathAbundanceEstimator::removeNoiseAndRenormalizeAbundances(PathClusterEstimates * path_cluster_estimates) const {

    const double noise_read_count = path_cluster_estimates->abundances(0, path_cluster_estimates->abundances.cols() - 1) * path_cluster_estimates->total_read_count;
    assert(noise_read_count <= path_cluster_estimates->total_read_count);

    path_cluster_estimates->posteriors.conservativeResize(1, path_cluster_estimates->posteriors.cols() - 1);
    path_cluster_estimates->abundances.conservativeResize(1, path_cluster_estimates->abundances.cols() - 1);

    assert(path_cluster_estimates->posteriors.cols() == path_cluster_estimates->paths.size());
    assert(path_cluster_estimates->abundances.cols() == path_cluster_estimates->paths.size());

    const double abundances_sum = path_cluster_estimates->abundances.sum();
    
    if (abundances_sum > 0) {

        path_cluster_estimates->abundances = path_cluster_estimates->abundances / abundances_sum;
    } 

    path_cluster_estimates->total_read_count -= noise_read_count;

    for (auto & abundance_samples: path_cluster_estimates->gibbs_read_count_samples) {

        assert(abundance_samples.path_ids.size() <= path_cluster_estimates->paths.size() + 1);
        assert(abundance_samples.samples.size() <= path_cluster_estimates->paths.size() + 1);
        
        assert(abundance_samples.path_ids.back() == path_cluster_estimates->paths.size());

        abundance_samples.path_ids.pop_back();
        abundance_samples.samples.pop_back();
    }
}

void PathAbundanceEstimator::updateEstimates(PathClusterEstimates * path_cluster_estimates, const PathClusterEstimates & new_path_cluster_estimates, const vector<uint32_t> & path_indices, const uint32_t sample_count) const {

    assert(new_path_cluster_estimates.posteriors.cols() == path_indices.size() + 1);
    assert(new_path_cluster_estimates.abundances.cols() == path_indices.size() + 1);

    assert(path_cluster_estimates->total_read_count == new_path_cluster_estimates.total_read_count);

   for (size_t i = 0; i < path_indices.size(); ++i) {

        path_cluster_estimates->posteriors(0, path_indices.at(i)) += (new_path_cluster_estimates.posteriors(0, i) * sample_count);            
        path_cluster_estimates->abundances(0, path_indices.at(i)) += (new_path_cluster_estimates.abundances(0, i) * sample_count);
    }

    path_cluster_estimates->posteriors(0, path_cluster_estimates->posteriors.cols() - 1) += (new_path_cluster_estimates.posteriors(0, path_indices.size()) * sample_count);            
    path_cluster_estimates->abundances(0, path_cluster_estimates->abundances.cols() - 1) += (new_path_cluster_estimates.abundances(0, path_indices.size()) * sample_count);  

    if (!new_path_cluster_estimates.gibbs_read_count_samples.empty()) {

       assert(new_path_cluster_estimates.gibbs_read_count_samples.size() == 1);
       path_cluster_estimates->gibbs_read_count_samples.emplace_back(move(new_path_cluster_estimates.gibbs_read_count_samples.front()));
    } 
}


MinimumPathAbundanceEstimator::MinimumPathAbundanceEstimator(const uint32_t max_em_its, const double min_em_conv, const uint32_t num_gibbs_samples, const uint32_t gibbs_thin_its, const double prob_precision) : PathAbundanceEstimator(max_em_its, min_em_conv, num_gibbs_samples, gibbs_thin_its, prob_precision) {}

void MinimumPathAbundanceEstimator::estimate(PathClusterEstimates * path_cluster_estimates, const vector<ReadPathProbabilities> & cluster_probs, mt19937 * mt_rng) {

    if (!cluster_probs.empty()) {

        Eigen::ColMatrixXd read_path_probs;
        Eigen::ColVectorXd noise_probs;
        Eigen::RowVectorXui read_counts;

        vector<uint32_t> path_ids(path_cluster_estimates->paths.size());
        iota(path_ids.begin(), path_ids.end(), 0);

        constructProbabilityMatrix(&read_path_probs, &noise_probs, &read_counts, cluster_probs, path_ids);      

        Eigen::ColMatrixXb read_path_cover = Eigen::ColMatrixXb::Zero(read_path_probs.rows(), read_path_probs.cols());
        Eigen::RowVectorXd path_weights = Eigen::RowVectorXd::Zero(read_path_probs.cols());

        for (size_t i = 0; i < read_path_probs.rows(); ++i) {

            if (doubleCompare(noise_probs(i), 1)) {

                read_counts(i) = 0;
            }

            for (auto & prob: cluster_probs.at(i).probabilities()) {

                assert(prob.second > 0);

                read_path_cover(i, prob.first) = true;
                path_weights(prob.first) += log(prob.second) * read_counts(i);                 
            }
        }

        path_weights *= -1;
        vector<uint32_t> min_path_cover = weightedMinimumPathCover(read_path_cover, read_counts, path_weights);

        if (!min_path_cover.empty()) {

            Eigen::ColMatrixXd min_path_read_path_probs;
            Eigen::ColVectorXd min_path_noise_probs;
            Eigen::RowVectorXui min_path_read_counts;

            constructProbabilityMatrix(&min_path_read_path_probs, &min_path_noise_probs, &min_path_read_counts, cluster_probs, min_path_cover);

            addNoiseAndNormalizeProbabilityMatrix(&min_path_read_path_probs, min_path_noise_probs);
            assert(min_path_read_path_probs.cols() >= 2);

            readCollapseProbabilityMatrix(&min_path_read_path_probs, &min_path_read_counts);

            PathClusterEstimates min_path_cluster_estimates;
            min_path_cluster_estimates.initEstimates(min_path_read_path_probs.cols(), 0, false);
            min_path_cluster_estimates.total_read_count = min_path_read_counts.sum();

            EMAbundanceEstimator(&min_path_cluster_estimates, min_path_read_path_probs, min_path_read_counts);
            assert(min_path_cluster_estimates.abundances.cols() == min_path_cover.size() + 1);            

            path_cluster_estimates->initEstimates(path_cluster_estimates->paths.size() + 1, 0, true);
            path_cluster_estimates->total_read_count = read_counts.sum();

            if (num_gibbs_samples > 0) {

                vector<CountSamples> * gibbs_read_count_samples = &(min_path_cluster_estimates.gibbs_read_count_samples);
                gibbs_read_count_samples->emplace_back(CountSamples());

                gibbs_read_count_samples->back().path_ids = min_path_cover;
                gibbs_read_count_samples->back().path_ids.emplace_back(path_cluster_estimates->abundances.cols() - 1);
                
                gibbs_read_count_samples->back().samples = vector<vector<double> >(min_path_cluster_estimates.abundances.cols(), vector<double>());

                gibbsReadCountSampler(&min_path_cluster_estimates, min_path_read_path_probs, min_path_read_counts, abundance_gibbs_gamma, mt_rng);
            }

            updateEstimates(path_cluster_estimates, min_path_cluster_estimates, min_path_cover, 1);    
            removeNoiseAndRenormalizeAbundances(path_cluster_estimates);

        } else {

            path_cluster_estimates->initEstimates(path_cluster_estimates->paths.size(), 0, true);
        }

    } else {

        path_cluster_estimates->initEstimates(path_cluster_estimates->paths.size(), 0, true);
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

        double max_weighted_read_path_cover = 0;
        int32_t max_weighted_read_path_cover_idx = -1;

        for (size_t i = 0; i < weighted_read_path_cover.size(); ++i) {

            if (weighted_read_path_cover(i) > max_weighted_read_path_cover) {

                max_weighted_read_path_cover = weighted_read_path_cover(i);
                max_weighted_read_path_cover_idx = i;
            }
        }

        assert(max_weighted_read_path_cover > 0);
        assert(max_weighted_read_path_cover_idx >= 0);

        min_path_cover.emplace_back(max_weighted_read_path_cover_idx);
        uncovered_read_counts = (uncovered_read_counts.array() * (!read_path_cover.col(max_weighted_read_path_cover_idx).transpose().array()).cast<uint32_t>()).matrix();
    }

    assert(min_path_cover.size() <= read_path_cover.cols());
    sort(min_path_cover.begin(), min_path_cover.end());

    return min_path_cover;
}


NestedPathAbundanceEstimator::NestedPathAbundanceEstimator(const uint32_t ploidy_in, const bool use_hap_gibbs_in, const uint32_t num_nested_samples_in, const uint32_t max_em_its, const double min_em_conv, const uint32_t num_gibbs_samples, const uint32_t gibbs_thin_its, const double prob_precision) : ploidy(ploidy_in), use_hap_gibbs(use_hap_gibbs_in), num_nested_samples(num_nested_samples_in), PathAbundanceEstimator(max_em_its, min_em_conv, num_gibbs_samples, gibbs_thin_its, prob_precision) {}

void NestedPathAbundanceEstimator::estimate(PathClusterEstimates * path_cluster_estimates, const vector<ReadPathProbabilities> & cluster_probs, mt19937 * mt_rng) {

    if (!cluster_probs.empty()) {

        auto path_groups = findPathOriginGroups(path_cluster_estimates->paths);

        vector<vector<uint32_t> > ploidy_path_indices_samples(num_nested_samples);

        for (auto & path_indices_samples: ploidy_path_indices_samples) {

            path_indices_samples.reserve(path_groups.size() * ploidy);
        }

        for (auto & group: path_groups) {    

            Eigen::ColMatrixXd group_read_path_probs;
            Eigen::ColVectorXd group_noise_probs;
            Eigen::RowVectorXui group_read_counts;        

            constructProbabilityMatrix(&group_read_path_probs, &group_noise_probs, &group_read_counts, cluster_probs, group);

            group_read_path_probs.conservativeResize(group_read_path_probs.rows(), group_read_path_probs.cols() + 1);
            group_read_path_probs.col(group_read_path_probs.cols() - 1) = group_noise_probs;

            readCollapseProbabilityMatrix(&group_read_path_probs, &group_read_counts);

            group_noise_probs = group_read_path_probs.col(group_read_path_probs.cols() - 1);
            group_read_path_probs.conservativeResize(group_read_path_probs.rows(), group_read_path_probs.cols() - 1);

            vector<uint32_t> group_path_counts;
            group_path_counts.reserve(group.size());

            for (size_t i = 0; i < group.size(); ++i) {

                group_path_counts.emplace_back(path_cluster_estimates->paths.at(group.at(i)).count);
            }

            PathClusterEstimates group_path_cluster_estimates;

            if (use_hap_gibbs) {

                estimatePathGroupPosteriorsGibbs(&group_path_cluster_estimates, group_read_path_probs, group_noise_probs, group_read_counts, group_path_counts, ploidy, mt_rng);

            } else {

                if (ploidy == 2) {

                    calculatePathGroupPosteriorsBounded(&group_path_cluster_estimates, group_read_path_probs, group_noise_probs, group_read_counts, group_path_counts, ploidy);

                } else {

                    calculatePathGroupPosteriorsFull(&group_path_cluster_estimates, group_read_path_probs, group_noise_probs, group_read_counts, group_path_counts, ploidy);                    
                }
            }

            samplePloidyPathIndices(&ploidy_path_indices_samples, group_path_cluster_estimates, group, mt_rng);
        }

        unordered_map<vector<uint32_t>, uint32_t> collapsed_ploidy_path_indices_samples;

        for (auto & path_samples: ploidy_path_indices_samples) {

            sort(path_samples.begin(), path_samples.end());

            auto collapsed_ploidy_path_indices_samples_it = collapsed_ploidy_path_indices_samples.emplace(path_samples, 0);
            collapsed_ploidy_path_indices_samples_it.first->second++;
        }

        path_cluster_estimates->initEstimates(path_cluster_estimates->paths.size() + 1, 0, true);

        for (auto & cluster_prob: cluster_probs) {

            path_cluster_estimates->total_read_count += cluster_prob.readCount();
        }

        for (auto & path_indices_sample: collapsed_ploidy_path_indices_samples) {

            assert(path_indices_sample.second > 0);

            Eigen::ColMatrixXd ploidy_read_path_probs;
            Eigen::ColVectorXd ploidy_noise_probs;
            Eigen::RowVectorXui ploidy_read_counts;

            constructProbabilityMatrix(&ploidy_read_path_probs, &ploidy_noise_probs, &ploidy_read_counts, cluster_probs, path_indices_sample.first);

            addNoiseAndNormalizeProbabilityMatrix(&ploidy_read_path_probs, ploidy_noise_probs);
            assert(ploidy_read_path_probs.cols() >= 2);

            readCollapseProbabilityMatrix(&ploidy_read_path_probs, &ploidy_read_counts);

            PathClusterEstimates ploidy_path_cluster_estimates;
            ploidy_path_cluster_estimates.initEstimates(ploidy_read_path_probs.cols(), 0, false);
            ploidy_path_cluster_estimates.total_read_count = ploidy_read_counts.sum();
           
            EMAbundanceEstimator(&ploidy_path_cluster_estimates, ploidy_read_path_probs, ploidy_read_counts);
            assert(ploidy_path_cluster_estimates.abundances.cols() == path_indices_sample.first.size() + 1);            

            if (num_gibbs_samples > 0) {

                vector<CountSamples> * gibbs_read_count_samples = &(ploidy_path_cluster_estimates.gibbs_read_count_samples);
                gibbs_read_count_samples->emplace_back(CountSamples());

                gibbs_read_count_samples->back().path_ids = path_indices_sample.first;
                gibbs_read_count_samples->back().path_ids.emplace_back(path_cluster_estimates->abundances.cols() - 1);
                
                gibbs_read_count_samples->back().samples = vector<vector<double> >(ploidy_path_cluster_estimates.abundances.cols(), vector<double>());

                for (uint32_t i = 0; i < path_indices_sample.second; ++i) {

                    gibbsReadCountSampler(&ploidy_path_cluster_estimates, ploidy_read_path_probs, ploidy_read_counts, abundance_gibbs_gamma, mt_rng);
                }
            }

            updateEstimates(path_cluster_estimates, ploidy_path_cluster_estimates, path_indices_sample.first, path_indices_sample.second);
        }

        for (size_t i = 0; i < path_cluster_estimates->abundances.cols(); ++i) {

            if (path_cluster_estimates->posteriors(0, i) > 0) {

                path_cluster_estimates->abundances(0, i) /= path_cluster_estimates->posteriors(0, i);
            }

            path_cluster_estimates->posteriors(0, i) /= num_nested_samples;
        }

        removeNoiseAndRenormalizeAbundances(path_cluster_estimates);

    } else {

        path_cluster_estimates->initEstimates(path_cluster_estimates->paths.size(), 0, true);
    }
}

vector<vector<uint32_t> > NestedPathAbundanceEstimator::findPathOriginGroups(const vector<PathInfo> & paths) const {

    vector<vector<uint32_t> > path_groups;
    unordered_map<string, uint32_t> path_group_indexes;

    for (size_t i = 0; i < paths.size(); ++i) {

        assert(paths.at(i).origin != "");

        auto path_group_indexes_it = path_group_indexes.emplace(paths.at(i).origin, path_group_indexes.size());

        if (path_group_indexes_it.second) {

            path_groups.emplace_back(vector<uint32_t>());
        }

        path_groups.at(path_group_indexes_it.first->second).emplace_back(i);
    }

    return path_groups;
}

void NestedPathAbundanceEstimator::samplePloidyPathIndices(vector<vector<uint32_t> > * ploidy_path_indices_samples, const PathClusterEstimates & group_path_cluster_estimates, const vector<uint32_t> & group, mt19937 * mt_rng) {

    const Eigen::RowVectorXd & posteriors = group_path_cluster_estimates.posteriors;
    assert(posteriors.cols() == group_path_cluster_estimates.path_groups.size());

    discrete_distribution<uint32_t> group_ploidy_path_sampler(posteriors.row(0).data(), posteriors.row(0).data() + posteriors.row(0).size());

    for (size_t i = 0; i < num_nested_samples; ++i) {

        vector<uint32_t> sampled_path_indices = group_path_cluster_estimates.path_groups.at(group_ploidy_path_sampler(*mt_rng));

        assert(!sampled_path_indices.empty());
        assert(sampled_path_indices.size() == ploidy);

        sort(sampled_path_indices.begin(), sampled_path_indices.end());

        ploidy_path_indices_samples->at(i).emplace_back(group.at(sampled_path_indices.front()));

        for (size_t j = 1; j < sampled_path_indices.size(); ++j) {

            if (ploidy_path_indices_samples->at(i).back() != group.at(sampled_path_indices.at(j))) {

                ploidy_path_indices_samples->at(i).emplace_back(group.at(sampled_path_indices.at(j)));
            }
        }
    }
}

