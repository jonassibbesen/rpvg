
#include <limits>
#include <chrono>

#include "sparsepp/spp.h"

#include "path_abundance_estimator.hpp"


const uint32_t min_em_conv_its = 10;
const double min_em_abundance = 1e-8;

const double abundance_gibbs_gamma = 1;
const double min_gibbs_abundance = 1e-8;

PathAbundanceEstimator::PathAbundanceEstimator(const uint32_t max_em_its_in, const double max_rel_em_conv_in, const uint32_t num_gibbs_samples_in, const uint32_t gibbs_thin_its_in, const double prob_precision) : max_em_its(max_em_its_in), max_rel_em_conv(max_rel_em_conv_in), num_gibbs_samples(num_gibbs_samples_in), gibbs_thin_its(gibbs_thin_its_in), PathEstimator(prob_precision) {}

void PathAbundanceEstimator::estimate(PathClusterEstimates * path_cluster_estimates, const vector<ReadPathProbabilities> & cluster_probs, mt19937 * mt_rng) {

    double debug_time1_0 = gbwt::readTimer();

    path_cluster_estimates->resetEstimates(path_cluster_estimates->paths.size(), 1);

    if (!cluster_probs.empty()) {

        Utils::ColMatrixXd read_path_probs;
        Utils::ColVectorXd noise_probs;
        Utils::RowVectorXd read_counts;

        double debug_time1_1 = gbwt::readTimer();

        if (path_cluster_estimates->out_debug) {

            #pragma omp critical
            {
                
                cerr << "DEBUG: INF1 " << omp_get_thread_num() << ": " << path_cluster_estimates->paths.size() << " " << cluster_probs.size() << " - " << gbwt::inGigabytes(gbwt::memoryUsage()) << " " << debug_time1_1 - debug_time1_0 << endl;
            }
        }

        constructProbabilityMatrix(&read_path_probs, &noise_probs, &read_counts, cluster_probs, path_cluster_estimates->paths.size());

        double debug_time1_2 = gbwt::readTimer();

        if (path_cluster_estimates->out_debug) {

            #pragma omp critical
            {
                
                cerr << "DEBUG: INF2 " << omp_get_thread_num() << ": " << path_cluster_estimates->paths.size() << " " << cluster_probs.size() << " " << read_path_probs.cols() << " " << read_path_probs.rows() << " - " << gbwt::inGigabytes(gbwt::memoryUsage()) << " " << debug_time1_2 - debug_time1_1 << endl;
            }
        }

        addNoiseAndNormalizeProbabilityMatrix(&read_path_probs, noise_probs);

        double debug_time1_3 = gbwt::readTimer();

        if (path_cluster_estimates->out_debug) {

            #pragma omp critical
            {
                
                cerr << "DEBUG: INF3 " << omp_get_thread_num() << ": " << path_cluster_estimates->paths.size() << " " << cluster_probs.size() << " " << read_path_probs.cols() << " " << read_path_probs.rows() << " - " << gbwt::inGigabytes(gbwt::memoryUsage()) << " " << debug_time1_3 - debug_time1_2 << endl;
            }
        }

        path_cluster_estimates->total_count = read_counts.sum();
        EMAbundanceEstimator(path_cluster_estimates, read_path_probs, read_counts);

        double debug_time1_4 = gbwt::readTimer();

        if (path_cluster_estimates->out_debug) {

            #pragma omp critical
            {
                
                cerr << "DEBUG: INF4 " << omp_get_thread_num() << ": " << path_cluster_estimates->paths.size() << " " << cluster_probs.size() << " " << read_path_probs.cols() << " " << read_path_probs.rows() << " " << path_cluster_estimates->path_group_sets.size() << " " << path_cluster_estimates->posteriors.size() << " " << path_cluster_estimates->abundances.size() << " - " << gbwt::inGigabytes(gbwt::memoryUsage()) << " " << debug_time1_4 - debug_time1_3 << endl;
            }
        }

        if (num_gibbs_samples > 0) {

            vector<CountSamples> * gibbs_read_count_samples = &(path_cluster_estimates->gibbs_read_count_samples);
            gibbs_read_count_samples->emplace_back(CountSamples());

            gibbs_read_count_samples->back().path_ids = vector<uint32_t>(path_cluster_estimates->path_group_sets.size());
            iota(gibbs_read_count_samples->back().path_ids.begin(), gibbs_read_count_samples->back().path_ids.end(), 0);

            gibbsReadCountSampler(path_cluster_estimates, read_path_probs, read_counts, abundance_gibbs_gamma, mt_rng, num_gibbs_samples);
        }
    } 
}

void PathAbundanceEstimator::EMAbundanceEstimator(PathClusterEstimates * path_cluster_estimates, const Utils::ColMatrixXd & read_path_probs, const Utils::RowVectorXd & read_counts) const {

    assert(!path_cluster_estimates->abundances.empty());

    assert(path_cluster_estimates->noise_count == 0);        
    assert(path_cluster_estimates->total_count > 0);

    Utils::RowVectorXd abundances = Eigen::RowVectorXd::Constant(1, path_cluster_estimates->abundances.size() + 1, 1 / static_cast<float>(path_cluster_estimates->abundances.size() + 1));
    Utils::RowVectorXd prev_abundances = abundances;

    uint32_t em_conv_its = 0;

    for (uint32_t i = 0; i < max_em_its; ++i) {

        Utils::ColMatrixXd read_posteriors = read_path_probs.array().rowwise() * abundances.array();
        read_posteriors = read_posteriors.array().colwise() / read_posteriors.rowwise().sum().array();

        abundances = read_counts * read_posteriors;
        abundances /= path_cluster_estimates->total_count;

        bool has_converged = true;

        for (size_t i = 0; i < abundances.cols(); ++i) {

            if (abundances(0, i) >= min_em_abundance) {

                auto rel_abundance_diff = fabs(abundances(0, i) - prev_abundances(0, i)) / abundances(0, i);

                if (rel_abundance_diff > max_rel_em_conv) {

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

        prev_abundances = abundances;
    }

    for (size_t i = 0; i < abundances.cols() - 1; ++i) {

        if (abundances(0, i) < min_em_abundance) {

            path_cluster_estimates->noise_count += abundances(0, i) * path_cluster_estimates->total_count;
            path_cluster_estimates->abundances.at(i) = 0;            

        } else {

            path_cluster_estimates->abundances.at(i) = abundances(0, i) * path_cluster_estimates->total_count;            
        }
    }

    path_cluster_estimates->noise_count += abundances(0, abundances.cols() - 1) * path_cluster_estimates->total_count;    
}

void PathAbundanceEstimator::EMAbundanceEstimator(PathClusterEstimates * path_cluster_estimates, const vector<ReadPathProbabilities> & cluster_probs) const {

    assert(!path_cluster_estimates->abundances.empty());

    assert(path_cluster_estimates->noise_count == 0);        
    assert(path_cluster_estimates->total_count > 0);

    vector<double> abundances(path_cluster_estimates->abundances.size() + 1, 1 / static_cast<float>(path_cluster_estimates->abundances.size() + 1));
    vector<double> prev_abundances = abundances;

    uint32_t em_conv_its = 0;

    for (uint32_t i = 0; i < max_em_its; ++i) {

        fill(abundances.begin(), abundances.end(), 0);

        for (auto & cluster_prob: cluster_probs) {

            double sum_prob = 0;
            assert(cluster_prob.isNoiseNorm());

            for (auto & path_probs: cluster_prob.pathProbs()) {

                for (auto & path: path_probs.second) {
        
                    assert(path < prev_abundances.size() - 1);
                    sum_prob += (prev_abundances.at(path) * path_probs.first);                
                }
            }

            sum_prob += prev_abundances.back() * cluster_prob.noiseProb();

            assert(!Utils::doubleCompare(sum_prob, 0));
            double norm_count = cluster_prob.readCount() / (sum_prob * path_cluster_estimates->total_count);

            for (auto & path_probs: cluster_prob.pathProbs()) {

                for (auto & path: path_probs.second) {
        
                    abundances.at(path) += (prev_abundances.at(path) * path_probs.first * norm_count);                
                }
            }   

            abundances.back() += (prev_abundances.back() * cluster_prob.noiseProb() * norm_count);       
        }

        bool has_converged = true;

        for (size_t i = 0; i < abundances.size(); ++i) {

            if (abundances.at(i) >= min_em_abundance) {

                auto rel_abundance_diff = fabs(abundances.at(i) - prev_abundances.at(i)) / abundances.at(i);

                if (rel_abundance_diff > max_rel_em_conv) {

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

        prev_abundances = abundances;
    }

    for (size_t i = 0; i < abundances.size() - 1; ++i) {

        if (abundances.at(i) < min_em_abundance) {

            path_cluster_estimates->noise_count += abundances.at(i) * path_cluster_estimates->total_count;
            path_cluster_estimates->abundances.at(i) = 0;            

        } else {

            path_cluster_estimates->abundances.at(i) = abundances.at(i) * path_cluster_estimates->total_count;            
        }
    }

    path_cluster_estimates->noise_count += abundances.back() * path_cluster_estimates->total_count;    
}

void PathAbundanceEstimator::gibbsReadCountSampler(PathClusterEstimates * path_cluster_estimates, const Utils::ColMatrixXd & read_path_probs, const Utils::RowVectorXd & read_counts, const double gamma, mt19937 * mt_rng, const uint32_t num_samples) const {

    assert(path_cluster_estimates->total_count > 0);

    assert(!path_cluster_estimates->gibbs_read_count_samples.empty());
    assert(path_cluster_estimates->gibbs_read_count_samples.back().path_ids.size() == path_cluster_estimates->abundances.size());

    assert(path_cluster_estimates->gibbs_read_count_samples.back().noise_samples.empty());
    path_cluster_estimates->gibbs_read_count_samples.back().noise_samples.reserve(num_samples);

    assert(path_cluster_estimates->gibbs_read_count_samples.back().abundance_samples.empty());
    path_cluster_estimates->gibbs_read_count_samples.back().abundance_samples.reserve(path_cluster_estimates->abundances.size() * num_samples);

    Utils::RowVectorXd gibbs_abundances = Eigen::RowVectorXd(1, path_cluster_estimates->abundances.size() + 1);

    for (size_t i = 0; i < gibbs_abundances.cols() - 1; ++i) {

        gibbs_abundances(0, i) = path_cluster_estimates->abundances.at(i) / path_cluster_estimates->total_count;
    }

    gibbs_abundances(0, gibbs_abundances.cols() - 1) = path_cluster_estimates->noise_count / path_cluster_estimates->total_count;    

    assert(Utils::doubleCompare(gibbs_abundances.sum(), 1));

    const uint32_t num_gibbs_its = num_samples * gibbs_thin_its;

    for (uint32_t gibbs_it = 1; gibbs_it <= num_gibbs_its; ++gibbs_it) {

        Utils::ColMatrixXd read_posteriors = read_path_probs.array().rowwise() * gibbs_abundances.array();
        read_posteriors = read_posteriors.array().colwise() / read_posteriors.rowwise().sum().array();

        vector<uint32_t> gibbs_path_read_counts(gibbs_abundances.cols(), 0);

        for (size_t i = 0; i < read_posteriors.rows(); ++i) {

            uint32_t row_reads_counts = read_counts(0, i);
            double row_sum_probs = 1;

            for (size_t j = 0; j < read_posteriors.cols(); ++j) {

                auto cur_prob = read_posteriors(i, j);

                if (cur_prob > 0) {

                    assert(row_sum_probs > 0);

                    binomial_distribution<uint32_t> path_read_count_sampler(row_reads_counts, min(1.0, cur_prob / row_sum_probs));
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

            path_cluster_estimates->gibbs_read_count_samples.back().noise_samples.emplace_back(0);

            for (size_t i = 0; i < gibbs_abundances.cols() - 1; ++i) {

                if (gibbs_abundances(0, i) < min_gibbs_abundance) {

                    path_cluster_estimates->gibbs_read_count_samples.back().noise_samples.back() += gibbs_abundances(0, i) * path_cluster_estimates->total_count;
                    path_cluster_estimates->gibbs_read_count_samples.back().abundance_samples.emplace_back(0);

                } else {

                    path_cluster_estimates->gibbs_read_count_samples.back().abundance_samples.emplace_back(gibbs_abundances(0, i) * path_cluster_estimates->total_count);
                }
            }

            path_cluster_estimates->gibbs_read_count_samples.back().noise_samples.back() += gibbs_abundances(0, gibbs_abundances.cols() - 1) * path_cluster_estimates->total_count;   
        }
    }
}


MinimumPathAbundanceEstimator::MinimumPathAbundanceEstimator(const uint32_t max_em_its, const double max_rel_em_conv, const uint32_t num_gibbs_samples, const uint32_t gibbs_thin_its, const double prob_precision) : PathAbundanceEstimator(max_em_its, max_rel_em_conv, num_gibbs_samples, gibbs_thin_its, prob_precision) {}

void MinimumPathAbundanceEstimator::estimate(PathClusterEstimates * path_cluster_estimates, const vector<ReadPathProbabilities> & cluster_probs, mt19937 * mt_rng) {

    path_cluster_estimates->resetEstimates(path_cluster_estimates->paths.size(), 1);

    if (!cluster_probs.empty()) {

        Utils::ColMatrixXb read_path_cover = Utils::ColMatrixXb::Zero(cluster_probs.size(), path_cluster_estimates->paths.size());
        Utils::RowVectorXd read_counts = Utils::RowVectorXd::Zero(cluster_probs.size());;
        Utils::RowVectorXd path_weights = Utils::RowVectorXd::Zero(path_cluster_estimates->paths.size());

        for (size_t i = 0; i < cluster_probs.size(); ++i) {
            
            assert(!cluster_probs.at(i).isNoiseNorm());
            auto noise_prob = cluster_probs.at(i).noiseProb();
    
            if (!Utils::doubleCompare(noise_prob, 1)) {

                read_counts(0, i) = cluster_probs.at(i).readCount();

                for (auto & path_probs: cluster_probs.at(i).pathProbs()) {

                    for (auto & path: path_probs.second) {

                        assert(path_probs.first > 0);

                        read_path_cover(i, path) = true;
                        path_weights(0, path) += log(path_probs.first * (1 - noise_prob)) * read_counts(0, i);  
                    }    
                }           
            }
        }

        path_weights *= -1;
        vector<uint32_t> min_path_cover = weightedMinimumPathCover(read_path_cover, read_counts, path_weights);

        if (!min_path_cover.empty()) {

            Utils::ColMatrixXd min_path_read_path_probs;
            Utils::ColVectorXd min_path_noise_probs;
            Utils::RowVectorXd min_path_read_counts;

            constructPartialProbabilityMatrix(&min_path_read_path_probs, &min_path_noise_probs, &min_path_read_counts, cluster_probs, min_path_cover, path_cluster_estimates->paths.size());
            
            PathClusterEstimates min_path_cluster_estimates;
            min_path_cluster_estimates.resetEstimates(min_path_read_path_probs.cols(), 1);

            addNoiseAndNormalizeProbabilityMatrix(&min_path_read_path_probs, min_path_noise_probs);
            readCollapseProbabilityMatrix(&min_path_read_path_probs, &min_path_read_counts);

            min_path_cluster_estimates.total_count = min_path_read_counts.sum();
            EMAbundanceEstimator(&min_path_cluster_estimates, min_path_read_path_probs, min_path_read_counts);

            if (num_gibbs_samples > 0) {

                vector<CountSamples> * gibbs_read_count_samples = &(min_path_cluster_estimates.gibbs_read_count_samples);
                gibbs_read_count_samples->emplace_back(CountSamples());

                gibbs_read_count_samples->back().path_ids = min_path_cover;                

                gibbsReadCountSampler(&min_path_cluster_estimates, min_path_read_path_probs, min_path_read_counts, abundance_gibbs_gamma, mt_rng, num_gibbs_samples);

                assert(min_path_cluster_estimates.gibbs_read_count_samples.size() == 1);
                path_cluster_estimates->gibbs_read_count_samples.emplace_back(move(min_path_cluster_estimates.gibbs_read_count_samples.front()));
            }

            assert(min_path_cluster_estimates.abundances.size() == min_path_cover.size());

            for (size_t i = 0; i < min_path_cover.size(); ++i) {

                path_cluster_estimates->abundances.at(min_path_cover.at(i)) += min_path_cluster_estimates.abundances.at(i);
            }

            path_cluster_estimates->noise_count = min_path_cluster_estimates.noise_count;
            path_cluster_estimates->total_count = min_path_cluster_estimates.total_count;
        } 
    } 
}

vector<uint32_t> MinimumPathAbundanceEstimator::weightedMinimumPathCover(const Utils::ColMatrixXb & read_path_cover, const Utils::RowVectorXd & read_counts, const Utils::RowVectorXd & path_weights) const {

    assert(read_path_cover.rows() == read_counts.cols());
    assert(read_path_cover.cols() == path_weights.cols());

    if (read_path_cover.cols() == 1) {

        return vector<uint32_t>({0});
    }

    auto uncovered_read_counts = read_counts;

    vector<uint32_t> min_path_cover;
    min_path_cover.reserve(read_path_cover.cols());

    while (uncovered_read_counts.maxCoeff() > 0) {

        Utils::RowVectorXd weighted_read_path_cover = (uncovered_read_counts * read_path_cover.cast<double>()).array() / path_weights.array();
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
        uncovered_read_counts = (uncovered_read_counts.array() * (!read_path_cover.col(max_weighted_read_path_cover_idx).transpose().array()).cast<double>()).matrix();
    }

    assert(min_path_cover.size() <= read_path_cover.cols());
    sort(min_path_cover.begin(), min_path_cover.end());

    return min_path_cover;
}

NestedPathAbundanceEstimator::NestedPathAbundanceEstimator(const uint32_t group_size_in, const double min_hap_prob_in, const bool infer_collapsed_in, const bool use_group_post_gibbs_in, const uint32_t max_em_its, const double max_rel_em_conv, const uint32_t num_gibbs_samples, const uint32_t gibbs_thin_its, const double prob_precision) : group_size(group_size_in), min_hap_prob(min_hap_prob_in), infer_collapsed(infer_collapsed_in), use_group_post_gibbs(use_group_post_gibbs_in), PathAbundanceEstimator(max_em_its, max_rel_em_conv, num_gibbs_samples, gibbs_thin_its, prob_precision) {}

void NestedPathAbundanceEstimator::estimate(PathClusterEstimates * path_cluster_estimates, const vector<ReadPathProbabilities> & cluster_probs, mt19937 * mt_rng) {

    if (infer_collapsed) {

        inferAbundancesCollapsedGroups(path_cluster_estimates, cluster_probs, mt_rng);
    
    } else {

        inferAbundancesIndependentGroups(path_cluster_estimates, cluster_probs, mt_rng);
    }
}

void NestedPathAbundanceEstimator::inferAbundancesIndependentGroups(PathClusterEstimates * path_cluster_estimates, const vector<ReadPathProbabilities> & cluster_probs, mt19937 * mt_rng) const {

    path_cluster_estimates->resetEstimates(0, 0);

    if (!cluster_probs.empty()) {

        auto path_groups = findPathGroups(path_cluster_estimates->paths);

        vector<PathClusterEstimates> group_path_cluster_estimates;
        group_path_cluster_estimates.reserve(path_groups.size());

        for (auto & group: path_groups) {    

            Utils::ColMatrixXd group_read_path_probs;
            Utils::ColVectorXd group_noise_probs;
            Utils::RowVectorXd group_read_counts;        

            constructPartialProbabilityMatrix(&group_read_path_probs, &group_noise_probs, &group_read_counts, cluster_probs, group, path_cluster_estimates->paths.size());

            addNoiseAndNormalizeProbabilityMatrix(&group_read_path_probs, group_noise_probs);
            readCollapseProbabilityMatrix(&group_read_path_probs, &group_read_counts);

            group_noise_probs = group_read_path_probs.col(group_read_path_probs.cols() - 1);
            group_read_path_probs.conservativeResize(group_read_path_probs.rows(), group_read_path_probs.cols() - 1);

            vector<uint32_t> group_path_counts;
            group_path_counts.reserve(group.size());

            for (size_t i = 0; i < group.size(); ++i) {

                group_path_counts.emplace_back(path_cluster_estimates->paths.at(group.at(i)).source_count);
            }

            group_path_cluster_estimates.emplace_back(PathClusterEstimates());

            if (use_group_post_gibbs) {

                estimatePathGroupPosteriorsGibbs(&(group_path_cluster_estimates.back()), group_read_path_probs, group_noise_probs, group_read_counts, group_path_counts, group_size, mt_rng);

            } else {

                if (group_size == 2) {

                    calculatePathGroupPosteriorsBounded(&(group_path_cluster_estimates.back()), group_read_path_probs, group_noise_probs, group_read_counts, group_path_counts, group_size, min_hap_prob);

                } else {

                    calculatePathGroupPosteriorsFull(&(group_path_cluster_estimates.back()), group_read_path_probs, group_noise_probs, group_read_counts, group_path_counts, group_size);                    
                }
            }
        }

        auto path_subset_samples = samplePathSubsetIndices(group_path_cluster_estimates, path_groups, mt_rng);
        inferPathSubsetAbundance(path_cluster_estimates, cluster_probs, mt_rng, path_subset_samples);
    } 
}

void NestedPathAbundanceEstimator::inferAbundancesCollapsedGroups(PathClusterEstimates * path_cluster_estimates, const vector<ReadPathProbabilities> & cluster_probs, mt19937 * mt_rng) const {

    path_cluster_estimates->resetEstimates(0, 0);

    if (!cluster_probs.empty()) {

        double debug_time0 = gbwt::readTimer();

        auto path_clusters = findPathClusters(path_cluster_estimates->paths);

        double debug_time1 = gbwt::readTimer();

        if (path_cluster_estimates->out_debug) {

            #pragma omp critical
            {
                
                cerr << "DEBUG: INF1 " << omp_get_thread_num() << ": " << path_cluster_estimates->paths.size() << " " << path_clusters.size() << " " << cluster_probs.size() << " - " << gbwt::inGigabytes(gbwt::memoryUsage()) << " " << debug_time1 - debug_time0 << endl;
            }
        }

        vector<pair<vector<vector<uint32_t> >, vector<uint32_t> > > path_source_groups;
        path_source_groups.reserve(path_clusters.size());

        vector<PathClusterEstimates> group_path_cluster_estimates;
        group_path_cluster_estimates.reserve(path_clusters.size());

        for (auto & cluster: path_clusters) {    

            // double debug_time1_0 = gbwt::readTimer();

            path_source_groups.emplace_back(findPathSourceGroups(path_cluster_estimates->paths, cluster));

            Utils::ColMatrixXd group_read_path_probs;
            Utils::ColVectorXd group_noise_probs;
            Utils::RowVectorXd group_read_counts;  

            // double debug_time1_1 = gbwt::readTimer();

            // if (path_cluster_estimates->out_debug) {

            //     #pragma omp critical
            //     {
                    
            //         cerr << "DEBUG: INF1.1 " << omp_get_thread_num() << ": " << path_cluster_estimates->paths.size() << " " << path_clusters.size() << " " << cluster_probs.size() << " " << cluster.size() << " " << path_source_groups.back().first.size() << " - " << gbwt::inGigabytes(gbwt::memoryUsage()) << " " << debug_time1_1 - debug_time1_0 << endl;
            //     }
            // }

            constructGroupedProbabilityMatrix(&group_read_path_probs, &group_noise_probs, &group_read_counts, cluster_probs, path_source_groups.back().first, path_cluster_estimates->paths.size());

            // double debug_time1_2 = gbwt::readTimer();

            // if (path_cluster_estimates->out_debug) {

            //     #pragma omp critical
            //     {
                    
            //         cerr << "DEBUG: INF1.2 " << omp_get_thread_num() << ": " << path_cluster_estimates->paths.size() << " " << path_clusters.size() << " " << cluster_probs.size() << " " << cluster.size() << " " << path_source_groups.back().first.size() << " " << group_read_path_probs.cols() << " " << group_read_path_probs.rows() << " - " << gbwt::inGigabytes(gbwt::memoryUsage()) << " " << debug_time1_2 - debug_time1_1 << endl;
            //     }
            // }
            
            addNoiseAndNormalizeProbabilityMatrix(&group_read_path_probs, group_noise_probs);
            readCollapseProbabilityMatrix(&group_read_path_probs, &group_read_counts);
            
            group_noise_probs = group_read_path_probs.col(group_read_path_probs.cols() - 1);
            group_read_path_probs.conservativeResize(group_read_path_probs.rows(), group_read_path_probs.cols() - 1);

            // double debug_time1_3 = gbwt::readTimer();

            // if (path_cluster_estimates->out_debug) {

            //     #pragma omp critical
            //     {
                    
            //         cerr << "DEBUG: INF1.3 " << omp_get_thread_num() << ": " << path_cluster_estimates->paths.size() << " " << path_clusters.size() << " " << cluster_probs.size() << " " << cluster.size() << " " << path_source_groups.back().first.size() << " " << group_read_path_probs.cols() << " " << group_read_path_probs.rows() << " - " << gbwt::inGigabytes(gbwt::memoryUsage()) << " " << debug_time1_3 - debug_time1_2 << endl;
            //     }
            // }

            group_path_cluster_estimates.emplace_back(PathClusterEstimates());

            if (use_group_post_gibbs) {

                estimatePathGroupPosteriorsGibbs(&(group_path_cluster_estimates.back()), group_read_path_probs, group_noise_probs, group_read_counts, path_source_groups.back().second, group_size, mt_rng);

            } else {

                if (group_size == 2) {

                    calculatePathGroupPosteriorsBounded(&(group_path_cluster_estimates.back()), group_read_path_probs, group_noise_probs, group_read_counts, path_source_groups.back().second, group_size, min_hap_prob);

                } else {

                    calculatePathGroupPosteriorsFull(&(group_path_cluster_estimates.back()), group_read_path_probs, group_noise_probs, group_read_counts, path_source_groups.back().second, group_size);                    
                }
            }

            // double debug_time1_4 = gbwt::readTimer();

            // if (path_cluster_estimates->out_debug) {

            //     #pragma omp critical
            //     {
                    
            //         cerr << "DEBUG: INF1.4 " << omp_get_thread_num() << ": " << path_cluster_estimates->paths.size() << " " << path_clusters.size() << " " << cluster_probs.size() << " " << cluster.size() << " " << path_source_groups.back().first.size() << " " << group_read_path_probs.cols() << " " << group_read_path_probs.rows() << " " << group_path_cluster_estimates.back().posteriors.size() << " - " << gbwt::inGigabytes(gbwt::memoryUsage()) << " " << debug_time1_4 - debug_time1_3 << endl;
            //     }
            // }
        }

        double debug_time2 = gbwt::readTimer();

        if (path_cluster_estimates->out_debug) {

            #pragma omp critical
            {
                
                cerr << "DEBUG: INF2 " << omp_get_thread_num() << ": " << path_cluster_estimates->paths.size() << " " << path_clusters.size() << " " << cluster_probs.size() << " - " << gbwt::inGigabytes(gbwt::memoryUsage()) << " " << debug_time2 - debug_time1 << endl;
            }
        }

        if (path_clusters.size() == 1) {
        
            auto path_subset_samples = selectPathSubsetIndices(group_path_cluster_estimates.front(), path_source_groups.front(), mt_rng);
            inferPathSubsetAbundance(path_cluster_estimates, cluster_probs, mt_rng, path_subset_samples);

        } else {

            auto path_subset_samples = samplePathSubsetIndices(group_path_cluster_estimates, path_source_groups, mt_rng);
        
            double debug_time3 = gbwt::readTimer();

            if (path_cluster_estimates->out_debug) {

                #pragma omp critical
                {
                    
                    cerr << "DEBUG: INF3 " << omp_get_thread_num() << ": " << path_cluster_estimates->paths.size() << " " << path_clusters.size() << " " << cluster_probs.size() << " " << path_subset_samples.size() << " - " << gbwt::inGigabytes(gbwt::memoryUsage()) << " " << debug_time3 - debug_time2 << endl;
                }
            }

            inferPathSubsetAbundance(path_cluster_estimates, cluster_probs, mt_rng, path_subset_samples);     

            if (path_cluster_estimates->out_debug) {

                #pragma omp critical
                {
                    
                    cerr << "DEBUG: INF4 " << omp_get_thread_num() << ": " << path_cluster_estimates->paths.size() << " " << path_clusters.size() << " " << cluster_probs.size() << " " << path_subset_samples.size() << " - " << gbwt::inGigabytes(gbwt::memoryUsage()) << " " << gbwt::readTimer() - debug_time3 << endl;
                }
            }
        }
    } 
}

vector<vector<uint32_t> > NestedPathAbundanceEstimator::findPathGroups(const vector<PathInfo> & paths) const {

    vector<vector<uint32_t> > path_groups;
    spp::sparse_hash_map<uint32_t, uint32_t> path_group_indexes;

    for (size_t i = 0; i < paths.size(); ++i) {

        auto path_group_indexes_it = path_group_indexes.emplace(paths.at(i).group_id, path_group_indexes.size());

        if (path_group_indexes_it.second) {

            path_groups.emplace_back(vector<uint32_t>());
        }

        path_groups.at(path_group_indexes_it.first->second).emplace_back(i);
    }

    return path_groups;
}

vector<vector<uint32_t> > NestedPathAbundanceEstimator::findPathClusters(const vector<PathInfo> & paths) const {

    vector<vector<uint32_t> > path_clusters;
    spp::sparse_hash_map<uint32_t, uint32_t> path_cluster_indexes;

    for (size_t i = 0; i < paths.size(); ++i) {

        auto path_cluster_indexes_it = path_cluster_indexes.emplace(paths.at(i).cluster_id, path_cluster_indexes.size());

        if (path_cluster_indexes_it.second) {

            path_clusters.emplace_back(vector<uint32_t>());
        }

        path_clusters.at(path_cluster_indexes_it.first->second).emplace_back(i);
    }

    return path_clusters;
}

pair<vector<vector<uint32_t> >, vector<uint32_t> > NestedPathAbundanceEstimator::findPathSourceGroups(const vector<PathInfo> & paths, const vector<uint32_t> & path_ids) const {

    spp::sparse_hash_map<uint32_t, vector<uint32_t> > source_id_paths;

    for (auto & path_id: path_ids) {

        for (auto & source_id: paths.at(path_id).source_ids) {

            auto source_id_paths_it = source_id_paths.emplace(source_id, vector<uint32_t>());
            source_id_paths_it.first->second.emplace_back(path_id);
        }
    }

    pair<vector<vector<uint32_t> >, vector<uint32_t> > path_source_groups;

    auto source_id_paths_it = source_id_paths.begin();

    while (source_id_paths_it != source_id_paths.end()) {

        if (source_id_paths_it->second.empty()) {

            ++source_id_paths_it;
            continue;
        }

        path_source_groups.second.emplace_back(1);

        auto source_id_paths_it2 = source_id_paths_it;
        ++source_id_paths_it2;

        while (source_id_paths_it2 != source_id_paths.end()) {

            if (!source_id_paths_it2->second.empty()) {

                if (source_id_paths_it->second == source_id_paths_it2->second) {

                    path_source_groups.second.back()++;
                    source_id_paths_it2->second.clear();
                } 
            }

            ++source_id_paths_it2;
        }

        assert(!source_id_paths_it->second.empty());
        path_source_groups.first.emplace_back(source_id_paths_it->second);

        source_id_paths_it->second.clear();
        ++source_id_paths_it;
    }

    assert(path_source_groups.first.size() == path_source_groups.second.size());
    return path_source_groups;
}

spp::sparse_hash_map<vector<uint32_t>, double> NestedPathAbundanceEstimator::samplePathSubsetIndices(const vector<PathClusterEstimates> & group_path_cluster_estimates, const vector<vector<uint32_t> > & path_groups, mt19937 * mt_rng) const {

    assert(group_path_cluster_estimates.size() == path_groups.size());

    vector<discrete_distribution<uint32_t> > path_group_set_samplers;
    path_group_set_samplers.reserve(group_path_cluster_estimates.size());

    for (auto & estimates: group_path_cluster_estimates) {

        assert(estimates.posteriors.size() == estimates.path_group_sets.size());
        path_group_set_samplers.emplace_back(estimates.posteriors.begin(), estimates.posteriors.end());
    }

    spp::sparse_hash_map<vector<uint32_t>, double> path_subset_samples;

    uint32_t num_samples = floor(1 / min_hap_prob);
    assert(num_samples > 0);

    for (size_t i = 0; i < num_samples; ++i) {

        vector<uint32_t> path_subset;

        for (size_t j = 0; j < group_path_cluster_estimates.size(); ++j) {

            vector<uint32_t> path_group_set = group_path_cluster_estimates.at(j).path_group_sets.at(path_group_set_samplers.at(j)(*mt_rng));

            assert(!path_group_set.empty());
            assert(path_group_set.size() == group_size);

            for (auto & group: path_group_set) {

                path_subset.emplace_back(path_groups.at(j).at(group));
            }
        }

        sort(path_subset.begin(), path_subset.end());

        auto path_subset_samples_it = path_subset_samples.emplace(path_subset, 0);
        path_subset_samples_it.first->second += 1 / static_cast<double>(num_samples);
    }

    return path_subset_samples;
}

spp::sparse_hash_map<vector<uint32_t>, double> NestedPathAbundanceEstimator::samplePathSubsetIndices(const vector<PathClusterEstimates> & group_path_cluster_estimates, const vector<pair<vector<vector<uint32_t> >, vector<uint32_t> > > & path_source_groups, mt19937 * mt_rng) const {

    assert(group_path_cluster_estimates.size() == path_source_groups.size());

    vector<discrete_distribution<uint32_t> > path_group_set_samplers;
    path_group_set_samplers.reserve(group_path_cluster_estimates.size());

    for (auto & estimates: group_path_cluster_estimates) {

        assert(estimates.posteriors.size() == estimates.path_group_sets.size());
        path_group_set_samplers.emplace_back(estimates.posteriors.begin(), estimates.posteriors.end());
    }

    spp::sparse_hash_map<vector<uint32_t>, double> path_subset_samples;

    uint32_t num_samples = floor(1 / min_hap_prob);
    assert(num_samples > 0);

    for (size_t i = 0; i < num_samples; ++i) {

        vector<uint32_t> path_subset;

        for (size_t j = 0; j < group_path_cluster_estimates.size(); ++j) {

            vector<uint32_t> path_group_set = group_path_cluster_estimates.at(j).path_group_sets.at(path_group_set_samplers.at(j)(*mt_rng));

            assert(!path_group_set.empty());
            assert(path_group_set.size() == group_size);

            for (auto & group: path_group_set) {

                path_subset.insert(path_subset.end(), path_source_groups.at(j).first.at(group).begin(), path_source_groups.at(j).first.at(group).end());
            }
        }

        sort(path_subset.begin(), path_subset.end());

        auto path_subset_samples_it = path_subset_samples.emplace(path_subset, 0);
        path_subset_samples_it.first->second += 1 / static_cast<double>(num_samples);
    }

    return path_subset_samples;
}

spp::sparse_hash_map<vector<uint32_t>, double> NestedPathAbundanceEstimator::selectPathSubsetIndices(const PathClusterEstimates & group_path_cluster_estimates, const pair<vector<vector<uint32_t> >, vector<uint32_t> > & path_source_groups, mt19937 * mt_rng) const {

    assert(group_path_cluster_estimates.posteriors.size() == group_path_cluster_estimates.path_group_sets.size());

    spp::sparse_hash_map<vector<uint32_t>, double> path_subset_samples;

    double sum_posterior = 0;
        
    for (size_t i = 0; i < group_path_cluster_estimates.posteriors.size(); ++i) {

        if (group_path_cluster_estimates.posteriors.at(i) >= min_hap_prob) {

            auto path_group_set = group_path_cluster_estimates.path_group_sets.at(i);

            assert(!path_group_set.empty());
            assert(path_group_set.size() == group_size);

            vector<uint32_t> path_subset;

            for (auto & group: path_group_set) {

                path_subset.insert(path_subset.end(), path_source_groups.first.at(group).begin(), path_source_groups.first.at(group).end());
            }

            sort(path_subset.begin(), path_subset.end());

            auto path_subset_samples_it = path_subset_samples.emplace(path_subset, 0);
            path_subset_samples_it.first->second += group_path_cluster_estimates.posteriors.at(i);

            sum_posterior += group_path_cluster_estimates.posteriors.at(i);
        }
    }

    for (auto & subset_sample: path_subset_samples) {

        subset_sample.second /= sum_posterior;
    }

    return path_subset_samples;
}

void NestedPathAbundanceEstimator::inferPathSubsetAbundance(PathClusterEstimates * path_cluster_estimates, const vector<ReadPathProbabilities> & cluster_probs, mt19937 * mt_rng, const spp::sparse_hash_map<vector<uint32_t>, double> & path_subset_samples) const {

    assert(path_cluster_estimates->noise_count == 0);
    assert(path_cluster_estimates->total_count == 0);

    for (auto & cluster_prob: cluster_probs) {

        path_cluster_estimates->total_count += cluster_prob.readCount();
    }

    spp::sparse_hash_map<vector<uint32_t>, pair<double, vector<double> > > path_group_estimates;

    double sum_hap_prob = 0;

    uint32_t subset_gibbs_samples = num_gibbs_samples;  
    double subset_gibbs_prob = 1;

    for (auto & path_subset: path_subset_samples) {

        double debug_time3_0 = gbwt::readTimer();

        if (path_subset.second < min_hap_prob) {

            continue;
        }
        
        sum_hap_prob += path_subset.second;

        assert(!path_subset.first.empty());
        assert(path_subset.second > 0);

        vector<uint32_t> collapsed_path_subset;
        collapsed_path_subset.reserve(path_subset.first.size());

        spp::sparse_hash_map<uint32_t, pair<uint32_t, uint32_t> > collapsed_path_subset_index;

        collapsed_path_subset.emplace_back(path_subset.first.front());
        collapsed_path_subset_index.emplace(path_subset.first.front(), make_pair(0, 1));

        for (size_t i = 1; i < path_subset.first.size(); ++i) {

            if (path_subset.first.at(i) != collapsed_path_subset.back()) {

                collapsed_path_subset.emplace_back(path_subset.first.at(i));
                collapsed_path_subset_index.emplace(path_subset.first.at(i), make_pair(collapsed_path_subset.size() - 1, 1));
            
            } else {

                collapsed_path_subset_index.at(path_subset.first.at(i)).second++;                
            }
        }


        double debug_time3_1 = gbwt::readTimer();

        if (path_cluster_estimates->out_debug) {

            #pragma omp critical
            {
                
                cerr << "DEBUG: INF3.1 " << omp_get_thread_num() << ": " << path_cluster_estimates->paths.size() << " " << path_subset_samples.size() << " " << cluster_probs.size() << " " << path_subset.first.size() << " " << path_subset.second << " " << collapsed_path_subset.size() << " - " << gbwt::inGigabytes(gbwt::memoryUsage()) << " " << debug_time3_1 - debug_time3_0 << endl;
            }
        }

        vector<ReadPathProbabilities> subset_cluster_probs;
        constructPartialClusterProbabilities(&subset_cluster_probs, cluster_probs, collapsed_path_subset, path_cluster_estimates->paths.size(), true);
        
        // Utils::ColMatrixXd subset_read_path_probs;
        // Utils::ColVectorXd subset_noise_probs;
        // Utils::RowVectorXd subset_read_counts;

        // constructPartialProbabilityMatrix(&subset_read_path_probs, &subset_noise_probs, &subset_read_counts, cluster_probs, collapsed_path_subset, path_cluster_estimates->paths.size());
     
        double debug_time3_2 = gbwt::readTimer();

        if (path_cluster_estimates->out_debug) {

            #pragma omp critical
            {
                
                cerr << "DEBUG: INF3.2 " << omp_get_thread_num() << ": " << path_cluster_estimates->paths.size() << " " << path_subset_samples.size() << " " << cluster_probs.size() << " " << path_subset.first.size() << " " << path_subset.second << " " << collapsed_path_subset.size() << " " << subset_cluster_probs.size() << " - " << gbwt::inGigabytes(gbwt::memoryUsage()) << " " << debug_time3_2 - debug_time3_1 << endl;
            }
        }

        // if (path_cluster_estimates->out_debug) {

        //     #pragma omp critical
        //     {
                
        //         cerr << "DEBUG: INF3.2 " << omp_get_thread_num() << ": " << path_cluster_estimates->paths.size() << " " << path_subset_samples.size() << " " << cluster_probs.size() << " " << path_subset.first.size() << " " << path_subset.second << " " << collapsed_path_subset.size() << " " << subset_read_path_probs.cols() << " " << subset_read_path_probs.rows() << " - " << gbwt::inGigabytes(gbwt::memoryUsage()) << " " << debug_time3_2 - debug_time3_1 << endl;
        //     }
        // }

        // PathClusterEstimates subset_path_cluster_estimates;
        // subset_path_cluster_estimates.resetEstimates(subset_read_path_probs.cols(), 1);
        
        // addNoiseAndNormalizeProbabilityMatrix(&subset_read_path_probs, subset_noise_probs);
        // readCollapseProbabilityMatrix(&subset_read_path_probs, &subset_read_counts);

        PathClusterEstimates subset_path_cluster_estimates;
        subset_path_cluster_estimates.resetEstimates(collapsed_path_subset.size(), 1);  

        double debug_time3_3 = gbwt::readTimer();

        if (path_cluster_estimates->out_debug) {

            #pragma omp critical
            {
                
                cerr << "DEBUG: INF3.3 " << omp_get_thread_num() << ": " << path_cluster_estimates->paths.size() << " " << path_subset_samples.size() << " " << cluster_probs.size() << " " << path_subset.first.size() << " " << path_subset.second << " " << collapsed_path_subset.size() << " " << subset_cluster_probs.size() << " - " << gbwt::inGigabytes(gbwt::memoryUsage()) << " " << debug_time3_3 - debug_time3_2 << endl;
            }
        }

        subset_path_cluster_estimates.total_count = 0;

        for (auto & cluster_prob: subset_cluster_probs) {

            subset_path_cluster_estimates.total_count += cluster_prob.readCount();
        }

        EMAbundanceEstimator(&subset_path_cluster_estimates, subset_cluster_probs);

        // if (path_cluster_estimates->out_debug) {

        //     #pragma omp critical
        //     {
                
        //         cerr << "DEBUG: INF3.3 " << omp_get_thread_num() << ": " << path_cluster_estimates->paths.size() << " " << path_subset_samples.size() << " " << cluster_probs.size() << " " << path_subset.first.size() << " " << path_subset.second << " " << collapsed_path_subset.size() << " " << subset_read_path_probs.cols() << " " << subset_read_path_probs.rows() << " - " << gbwt::inGigabytes(gbwt::memoryUsage()) << " " << debug_time3_3 - debug_time3_2 << endl;

        //         cerr << subset_read_path_probs << endl;
        //     }
        // }

        // subset_path_cluster_estimates.total_count = subset_read_counts.sum();
        // EMAbundanceEstimator(&subset_path_cluster_estimates, subset_read_path_probs, subset_read_counts);

        double debug_time3_4 = gbwt::readTimer();

        if (path_cluster_estimates->out_debug) {

            #pragma omp critical
            {
                
                cerr << "DEBUG: INF3.4 " << omp_get_thread_num() << ": " << path_cluster_estimates->paths.size() << " " << path_subset_samples.size() << " " << cluster_probs.size() << " " << path_subset.first.size() << " " << path_subset.second << " " << collapsed_path_subset.size() << " " << subset_cluster_probs.size() << " - " << gbwt::inGigabytes(gbwt::memoryUsage()) << " " << debug_time3_4 - debug_time3_3 << endl;
            }
        }
        
        // if (path_cluster_estimates->out_debug) {

        //     #pragma omp critical
        //     {
                
        //         cerr << "DEBUG: INF3.4 " << omp_get_thread_num() << ": " << path_cluster_estimates->paths.size() << " " << path_subset_samples.size() << " " << cluster_probs.size() << " " << path_subset.first.size() << " " << path_subset.second << " " << collapsed_path_subset.size() << " " << subset_read_path_probs.cols() << " " << subset_read_path_probs.rows() << " - " << gbwt::inGigabytes(gbwt::memoryUsage()) << " " << debug_time3_4 - debug_time3_3 << endl;
        //     }
        // }
        
        assert(subset_path_cluster_estimates.abundances.size() == collapsed_path_subset.size());            

        // if (subset_gibbs_samples > 0) {

        //     assert(subset_gibbs_prob > 0);

        //     binomial_distribution<uint32_t> path_read_count_sampler(subset_gibbs_samples, min(1.0, path_subset.second / subset_gibbs_prob));
        //     uint32_t cur_subset_gibbs_samples = path_read_count_sampler(*mt_rng);
            
        //     subset_gibbs_samples -= cur_subset_gibbs_samples;
        //     subset_gibbs_prob -= path_subset.second;

        //     if (cur_subset_gibbs_samples > 0) {

        //         vector<CountSamples> * gibbs_read_count_samples = &(subset_path_cluster_estimates.gibbs_read_count_samples);
        //         gibbs_read_count_samples->emplace_back(CountSamples());

        //         gibbs_read_count_samples->back().path_ids = collapsed_path_subset;            

        //         gibbsReadCountSampler(&subset_path_cluster_estimates, subset_read_path_probs, subset_read_counts, abundance_gibbs_gamma, mt_rng, cur_subset_gibbs_samples);

        //         assert(subset_path_cluster_estimates.gibbs_read_count_samples.size() == 1);
        //         path_cluster_estimates->gibbs_read_count_samples.emplace_back(move(subset_path_cluster_estimates.gibbs_read_count_samples.front()));
        //     }
        // }

        assert(subset_path_cluster_estimates.posteriors.size() == collapsed_path_subset.size());
        assert(subset_path_cluster_estimates.abundances.size() == collapsed_path_subset.size());

        path_cluster_estimates->noise_count += subset_path_cluster_estimates.noise_count * path_subset.second;
        assert(path_cluster_estimates->total_count == subset_path_cluster_estimates.total_count);

        spp::sparse_hash_map<uint32_t, vector<uint32_t> > subset_path_group_index;
        
        for (auto & path: path_subset.first) {

            auto subset_path_group_index_it = subset_path_group_index.emplace(path_cluster_estimates->paths.at(path).group_id, vector<uint32_t>());
            subset_path_group_index_it.first->second.emplace_back(path);
        }

        for (auto & path_group: subset_path_group_index) {

            assert(path_group.second.size() <= group_size);

            auto path_group_estimates_it = path_group_estimates.emplace(path_group.second, pair<uint32_t, vector<double> >(0, vector<double>(path_group.second.size(), 0)));
            path_group_estimates_it.first->second.first += path_subset.second;

            for (size_t i = 0; i < path_group.second.size(); ++i) {

                auto collapsed_path_subset_index_it = collapsed_path_subset_index.find(path_group.second.at(i));
                assert(collapsed_path_subset_index_it != collapsed_path_subset_index.end());

                path_group_estimates_it.first->second.second.at(i) += (subset_path_cluster_estimates.abundances.at(collapsed_path_subset_index_it->second.first) * path_subset.second / collapsed_path_subset_index_it->second.second);
            }          
        }

        if (path_cluster_estimates->out_debug) {

            #pragma omp critical
            {
                
                cerr << "DEBUG: INF3.5 " << omp_get_thread_num() << ": " << path_cluster_estimates->paths.size() << " " << path_subset_samples.size() << " " << cluster_probs.size() << " " << path_subset.first.size() << " " << path_subset.second << " " << collapsed_path_subset.size() << " " << subset_cluster_probs.size() << " - " << gbwt::inGigabytes(gbwt::memoryUsage()) << " " << gbwt::readTimer() - debug_time3_4 << endl;
            }
        }

        // if (path_cluster_estimates->out_debug) {

        //     #pragma omp critical
        //     {
                
        //         cerr << "DEBUG: INF3.5 " << omp_get_thread_num() << ": " << path_cluster_estimates->paths.size() << " " << path_subset_samples.size() << " " << cluster_probs.size() << " " << path_subset.first.size() << " " << path_subset.second << " " << collapsed_path_subset.size() << " " << subset_read_path_probs.cols() << " " << subset_read_path_probs.rows() << " - " << gbwt::inGigabytes(gbwt::memoryUsage()) << " " << gbwt::readTimer() - debug_time3_4 << endl;
        //     }
        // }
    }

    assert(path_cluster_estimates->path_group_sets.empty());
    assert(path_cluster_estimates->posteriors.empty());
    assert(path_cluster_estimates->abundances.empty());

    path_cluster_estimates->path_group_sets.reserve(path_group_estimates.size());
    path_cluster_estimates->posteriors.reserve(path_group_estimates.size());
    path_cluster_estimates->abundances.reserve(path_group_estimates.size() * group_size);

    for (auto & estimates: path_group_estimates) {

        assert(estimates.first.size() <= group_size);
        assert(estimates.first.size() == estimates.second.second.size());

        path_cluster_estimates->path_group_sets.emplace_back(estimates.first);
        path_cluster_estimates->posteriors.emplace_back(estimates.second.first);
        path_cluster_estimates->abundances.insert(path_cluster_estimates->abundances.end(), estimates.second.second.begin(), estimates.second.second.end());
    }

    assert(sum_hap_prob < 1 || Utils::doubleCompare(sum_hap_prob, 1));
    path_cluster_estimates->noise_count += (1 - sum_hap_prob) * path_cluster_estimates->total_count;
}

