
#include <limits>
#include <chrono>

#include "path_abundance_estimator.hpp"


const uint32_t min_em_conv_its = 10;
const double min_abundances = pow(10, -8);

PathAbundanceEstimator::PathAbundanceEstimator(const uint32_t max_em_its_in, const double min_em_conv, const double prob_precision) : max_em_its(max_em_its_in), em_conv_min_exp(min_em_conv), em_conv_max_rel_diff(min_em_conv), PathEstimator(prob_precision) {}

void PathAbundanceEstimator::estimate(PathClusterEstimates * path_cluster_estimates, const vector<ReadPathProbabilities> & cluster_probs) {

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

        EMAbundanceEstimator(path_cluster_estimates, read_path_probs, read_counts);
        removeNoiseAndRenormalizeAbundances(path_cluster_estimates);

    } else {

        path_cluster_estimates->initEstimates(path_cluster_estimates->paths.size(), 0, true);
    }
}

void PathAbundanceEstimator::EMAbundanceEstimator(PathClusterEstimates * path_cluster_estimates, const Eigen::ColMatrixXd & read_path_probs, const Eigen::RowVectorXui & read_counts) const {

    path_cluster_estimates->read_count = read_counts.sum();
    assert(path_cluster_estimates->read_count > 0);

    Eigen::RowVectorXd prev_abundances = path_cluster_estimates->abundances;
    uint32_t em_conv_its = 0;

    for (size_t i = 0; i < max_em_its; ++i) {

        Eigen::ColMatrixXd read_posteriors = read_path_probs.array().rowwise() * path_cluster_estimates->abundances.array();
        read_posteriors = read_posteriors.array().colwise() / read_posteriors.rowwise().sum().array();

        path_cluster_estimates->abundances = read_counts.cast<double>() * read_posteriors;
        path_cluster_estimates->abundances /= path_cluster_estimates->read_count;

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

        if (path_cluster_estimates->abundances(0, i) < min_abundances) {

            path_cluster_estimates->abundances(0, i) = 0;                    
        } 

        abundances_sum += path_cluster_estimates->abundances(0, i);
    }

    if (abundances_sum > 0) {

        path_cluster_estimates->abundances = path_cluster_estimates->abundances / abundances_sum;
    }
}

void PathAbundanceEstimator::removeNoiseAndRenormalizeAbundances(PathClusterEstimates * path_cluster_estimates) const {

    const double noise_read_count = path_cluster_estimates->abundances(0, path_cluster_estimates->abundances.cols() - 1) * path_cluster_estimates->read_count;

    if (!(noise_read_count <= path_cluster_estimates->read_count)) {

        cerr << noise_read_count << endl;
        cerr << path_cluster_estimates->read_count << endl;
        cerr << path_cluster_estimates->abundances << endl;
    }
    
    assert(noise_read_count <= path_cluster_estimates->read_count);

    path_cluster_estimates->posteriors.conservativeResize(1, path_cluster_estimates->posteriors.cols() - 1);
    path_cluster_estimates->abundances.conservativeResize(1, path_cluster_estimates->abundances.cols() - 1);

    const double abundances_sum = path_cluster_estimates->abundances.sum();  

    if (abundances_sum > 0) {

        path_cluster_estimates->abundances = path_cluster_estimates->abundances / abundances_sum;
    } 

    path_cluster_estimates->read_count -= noise_read_count;
}

MinimumPathAbundanceEstimator::MinimumPathAbundanceEstimator(const uint32_t max_em_its, const double min_em_conv, const double prob_precision) : PathAbundanceEstimator(max_em_its, min_em_conv, prob_precision) {}

void MinimumPathAbundanceEstimator::estimate(PathClusterEstimates * path_cluster_estimates, const vector<ReadPathProbabilities> & cluster_probs) {

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
            
            EMAbundanceEstimator(&min_path_cluster_estimates, min_path_read_path_probs, min_path_read_counts);

            path_cluster_estimates->initEstimates(path_cluster_estimates->paths.size() + 1, 0, true);
            path_cluster_estimates->read_count = min_path_read_counts.sum();

            for (size_t i = 0; i < min_path_cover.size(); i++) {

                path_cluster_estimates->posteriors(0, min_path_cover.at(i)) = min_path_cluster_estimates.posteriors(0, i);
                path_cluster_estimates->abundances(0, min_path_cover.at(i)) = min_path_cluster_estimates.abundances(0, i);
            }

            assert(min_path_cluster_estimates.posteriors.cols() == min_path_cover.size() + 1);

            path_cluster_estimates->posteriors(0, min_path_cover.size()) = min_path_cluster_estimates.posteriors(0, min_path_cover.size());
            path_cluster_estimates->abundances(0, min_path_cover.size()) = min_path_cluster_estimates.abundances(0, min_path_cover.size());  
                      
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

NestedPathAbundanceEstimator::NestedPathAbundanceEstimator(const uint32_t num_nested_its_in, const uint32_t ploidy_in, const bool use_exact_in, const uint32_t rng_seed, const uint32_t max_em_its, const double min_em_conv, const double prob_precision) : num_nested_its(num_nested_its_in), ploidy(ploidy_in), use_exact(use_exact_in), PathAbundanceEstimator(max_em_its, min_em_conv, prob_precision) {

    mt_rng = mt19937(rng_seed);
}

void NestedPathAbundanceEstimator::estimate(PathClusterEstimates * path_cluster_estimates, const vector<ReadPathProbabilities> & cluster_probs) {

    if (!cluster_probs.empty()) {

        auto path_groups = findPathOriginGroups(path_cluster_estimates->paths);

        vector<vector<uint32_t> > ploidy_path_indices_samples(num_nested_its);

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

            if (use_exact) {

                calculatePathGroupPosteriors(&group_path_cluster_estimates, group_read_path_probs, group_noise_probs, group_read_counts, group_path_counts, ploidy);

            } else {

                estimatePathGroupPosteriorsGibbs(&group_path_cluster_estimates, group_read_path_probs, group_noise_probs, group_read_counts, group_path_counts, ploidy, &mt_rng);
            }

            // Debug start

            if (path_cluster_estimates->paths.at(group.front()).origin == "ENST00000646664.1" || 
                path_cluster_estimates->paths.at(group.front()).origin == "ENST00000227378.7" || 
                path_cluster_estimates->paths.at(group.front()).origin == "ENST00000514057.1" || 
                path_cluster_estimates->paths.at(group.front()).origin == "ENST00000394667.7" || 
                path_cluster_estimates->paths.at(group.front()).origin == "ENST00000253788.11") {

                stringstream debug_stream;
                debug_stream << "\n######\n" << endl;

                for (size_t i = 0; i < group.size(); ++i) {

                    debug_stream << path_cluster_estimates->paths.at(group.at(i)).name << "\t" << path_cluster_estimates->paths.at(group.at(i)).count << endl;
                }

                debug_stream << endl;

                assert(group_path_cluster_estimates.path_groups.size() == group_path_cluster_estimates.posteriors.cols());

                for (size_t i = 0; i < group_path_cluster_estimates.path_groups.size(); ++i) {

                    if (group_path_cluster_estimates.posteriors(i) > 0) {

                        for (auto & id: group_path_cluster_estimates.path_groups.at(i)) {

                            debug_stream << path_cluster_estimates->paths.at(group.at(id)).name << ",";
                        }                    

                        debug_stream << "\t" << group_path_cluster_estimates.posteriors(i) << endl;
                    }
                }

                cerr << debug_stream.str() << endl;
            }

            // Debug end

            samplePloidyPathIndices(&ploidy_path_indices_samples, group_path_cluster_estimates, group);
        }

        unordered_map<vector<uint32_t>, uint32_t> collapsed_ploidy_path_indices_samples;

        for (auto & path_samples: ploidy_path_indices_samples) {

            sort(path_samples.begin(), path_samples.end());

            auto collapsed_ploidy_path_indices_samples_it = collapsed_ploidy_path_indices_samples.emplace(path_samples, 0);
            collapsed_ploidy_path_indices_samples_it.first->second++;
        }

        path_cluster_estimates->initEstimates(path_cluster_estimates->paths.size() + 1, 0, true);
        bool is_first = true;

        for (auto & path_indices_sample: collapsed_ploidy_path_indices_samples) {

            assert(path_indices_sample.second > 0);

            Eigen::ColMatrixXd ploidy_read_path_probs;
            Eigen::ColVectorXd ploidy_noise_probs;
            Eigen::RowVectorXui ploidy_read_counts;

            constructProbabilityMatrix(&ploidy_read_path_probs, &ploidy_noise_probs, &ploidy_read_counts, cluster_probs, path_indices_sample.first);

            addNoiseAndNormalizeProbabilityMatrix(&ploidy_read_path_probs, ploidy_noise_probs);
            assert(ploidy_read_path_probs.cols() >= 2);

            readCollapseProbabilityMatrix(&ploidy_read_path_probs, &ploidy_read_counts);

            if (is_first) {

                path_cluster_estimates->read_count = ploidy_read_counts.sum();
                is_first = false;
            } 

            PathClusterEstimates ploidy_path_cluster_estimates;
            ploidy_path_cluster_estimates.initEstimates(ploidy_read_path_probs.cols(), 0, false);
            
            EMAbundanceEstimator(&ploidy_path_cluster_estimates, ploidy_read_path_probs, ploidy_read_counts);
            updateEstimates(path_cluster_estimates, ploidy_path_cluster_estimates, path_indices_sample.first, path_indices_sample.second);
        }

        for (size_t i = 0; i < path_cluster_estimates->abundances.cols(); ++i) {

            if (path_cluster_estimates->posteriors(0, i) > 0) {

                path_cluster_estimates->abundances(0, i) /= path_cluster_estimates->posteriors(0, i);
            }

            path_cluster_estimates->posteriors(0, i) /= num_nested_its;
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

void NestedPathAbundanceEstimator::samplePloidyPathIndices(vector<vector<uint32_t> > * ploidy_path_indices_samples, const PathClusterEstimates & group_path_cluster_estimates, const vector<uint32_t> & group) {

    discrete_distribution<uint32_t> group_ploidy_path_sampler(group_path_cluster_estimates.posteriors.row(0).begin(), group_path_cluster_estimates.posteriors.row(0).end());

    assert(group_path_cluster_estimates.posteriors.cols() == group_path_cluster_estimates.path_groups.size());

    for (size_t i = 0; i < num_nested_its; ++i) {

        vector<uint32_t> sampled_path_indices = group_path_cluster_estimates.path_groups.at(group_ploidy_path_sampler(mt_rng));

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

void NestedPathAbundanceEstimator::updateEstimates(PathClusterEstimates * path_cluster_estimates, const PathClusterEstimates & new_path_cluster_estimates, const vector<uint32_t> & path_indices, const uint32_t sample_count) const {

   for (size_t i = 0; i < path_indices.size(); ++i) {

        if (new_path_cluster_estimates.posteriors(0, i) > 0) {

            assert(doubleCompare(new_path_cluster_estimates.posteriors(0, i), 1));

            path_cluster_estimates->posteriors(0, path_indices.at(i)) += (new_path_cluster_estimates.posteriors(0, i) * sample_count);            
            path_cluster_estimates->abundances(0, path_indices.at(i)) += (new_path_cluster_estimates.abundances(0, i) * sample_count);
        }
    }

    assert(new_path_cluster_estimates.posteriors.cols() == path_indices.size() + 1);

    if (new_path_cluster_estimates.posteriors(path_indices.size()) > 0) {

        path_cluster_estimates->posteriors(0, path_cluster_estimates->posteriors.cols() - 1) += (new_path_cluster_estimates.posteriors(0, path_indices.size()) * sample_count);
        path_cluster_estimates->abundances(0, path_cluster_estimates->abundances.cols() - 1) += (new_path_cluster_estimates.abundances(0, path_indices.size()) * sample_count);  
    }
}
