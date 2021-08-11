
#include <limits>
#include <chrono>

#include "sparsepp/spp.h"

#include "path_abundance_estimator.hpp"


const uint32_t min_em_conv_its = 10;
const double min_em_abundance = 1e-8;

const double abundance_gibbs_gamma = 1;

const uint32_t min_rel_likelihood_scaling = 1e4;

PathAbundanceEstimator::PathAbundanceEstimator(const uint32_t max_em_its_in, const double max_rel_em_conv_in, const uint32_t num_gibbs_samples_in, const uint32_t gibbs_thin_its_in, const double prob_precision) : max_em_its(max_em_its_in), max_rel_em_conv(max_rel_em_conv_in), num_gibbs_samples(num_gibbs_samples_in), gibbs_thin_its(gibbs_thin_its_in), PathEstimator(prob_precision) {}

void PathAbundanceEstimator::estimate(PathClusterEstimates * path_cluster_estimates, const vector<ReadPathProbabilities> & cluster_probs, mt19937 * mt_rng) {

    if (!cluster_probs.empty()) {

        Utils::ColMatrixXd read_path_probs;
        Utils::ColVectorXd noise_probs;
        Utils::RowVectorXd read_counts;

        constructProbabilityMatrix(&read_path_probs, &noise_probs, &read_counts, cluster_probs, path_cluster_estimates->paths.size());
        detractNoiseAndNormalizeProbabilityMatrix(&read_path_probs, &noise_probs, &read_counts);

        if (read_path_probs.rows() == 0) {

            assert(noise_probs.rows() == 0);
            assert(read_counts.cols() == 0);

            path_cluster_estimates->initEstimates(path_cluster_estimates->paths.size(), 0, true);
            return;
        }

        const double total_read_count = read_counts.sum();
        assert(total_read_count > 0);

        path_cluster_estimates->initEstimates(path_cluster_estimates->paths.size(), 0, false);
        EMAbundanceEstimator(path_cluster_estimates, read_path_probs, read_counts, total_read_count);

        if (num_gibbs_samples > 0) {

            vector<CountSamples> * gibbs_read_count_samples = &(path_cluster_estimates->gibbs_read_count_samples);
            gibbs_read_count_samples->emplace_back(CountSamples());

            gibbs_read_count_samples->back().path_ids = vector<uint32_t>(path_cluster_estimates->abundances.cols());
            iota(gibbs_read_count_samples->back().path_ids.begin(), gibbs_read_count_samples->back().path_ids.end(), 0);

            gibbs_read_count_samples->back().samples = vector<vector<double> >(path_cluster_estimates->abundances.cols(), vector<double>());

            gibbsReadCountSampler(path_cluster_estimates, read_path_probs, read_counts, total_read_count, abundance_gibbs_gamma, mt_rng);
        }

        path_cluster_estimates->abundances *= total_read_count;

    } else {

        path_cluster_estimates->initEstimates(path_cluster_estimates->paths.size(), 0, true);
    }
}

void PathAbundanceEstimator::EMAbundanceEstimator(PathClusterEstimates * path_cluster_estimates, const Utils::ColMatrixXd & read_path_probs, const Utils::RowVectorXd & read_counts, const double total_read_count) const {

    Utils::RowVectorXd prev_abundances = path_cluster_estimates->abundances;
    uint32_t em_conv_its = 0;

    for (uint32_t i = 0; i < max_em_its; ++i) {

        Utils::ColMatrixXd read_posteriors = read_path_probs.array().rowwise() * path_cluster_estimates->abundances.array();
        read_posteriors = read_posteriors.array().colwise() / read_posteriors.rowwise().sum().array();

        path_cluster_estimates->abundances = read_counts * read_posteriors;
        path_cluster_estimates->abundances /= total_read_count;

        bool has_converged = true;

        for (size_t i = 0; i < path_cluster_estimates->abundances.cols(); ++i) {

            if (path_cluster_estimates->abundances(0, i) >= min_em_abundance) {

                auto rel_abundance_diff = fabs(path_cluster_estimates->abundances(0, i) - prev_abundances(0, i)) / path_cluster_estimates->abundances(0, i);

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

        prev_abundances = path_cluster_estimates->abundances;
    }

    double abundances_sum = 0;

    for (size_t i = 0; i < path_cluster_estimates->abundances.cols(); ++i) {

        if (path_cluster_estimates->abundances(0, i) < min_em_abundance) {

            path_cluster_estimates->abundances(0, i) = 0;                    
        } 

        abundances_sum += path_cluster_estimates->abundances(0, i);
    }

    if (abundances_sum > 0) {

        path_cluster_estimates->abundances = path_cluster_estimates->abundances / abundances_sum;
    }
}

void PathAbundanceEstimator::gibbsReadCountSampler(PathClusterEstimates * path_cluster_estimates, const Utils::ColMatrixXd & read_path_probs, const Utils::RowVectorXd & read_counts, const double total_read_count, const double gamma, mt19937 * mt_rng) const {

    assert(!path_cluster_estimates->gibbs_read_count_samples.empty());
    assert(path_cluster_estimates->gibbs_read_count_samples.back().path_ids.size() == path_cluster_estimates->abundances.cols());
    assert(path_cluster_estimates->gibbs_read_count_samples.back().samples.size() == path_cluster_estimates->abundances.cols());

    assert(Utils::doubleCompare(path_cluster_estimates->abundances.sum(), 1));
    Utils::RowVectorXd gibbs_abundances = path_cluster_estimates->abundances;

    const uint32_t num_gibbs_its = num_gibbs_samples * gibbs_thin_its;

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

            for (size_t i = 0; i < gibbs_abundances.cols(); ++i) {

                path_cluster_estimates->gibbs_read_count_samples.back().samples.at(i).emplace_back(gibbs_abundances(0, i) * total_read_count);
            }
        }
    }
}

void PathAbundanceEstimator::updateEstimates(PathClusterEstimates * path_cluster_estimates, const PathClusterEstimates & new_path_cluster_estimates, const vector<uint32_t> & path_indices, const uint32_t sample_count) const {  

    assert(new_path_cluster_estimates.abundances.cols() == path_indices.size());

   for (size_t i = 0; i < path_indices.size(); ++i) {

        path_cluster_estimates->abundances(0, path_indices.at(i)) += (new_path_cluster_estimates.abundances(0, i) * sample_count);
    }

    if (!new_path_cluster_estimates.gibbs_read_count_samples.empty()) {

       assert(new_path_cluster_estimates.gibbs_read_count_samples.size() == 1);
       path_cluster_estimates->gibbs_read_count_samples.emplace_back(move(new_path_cluster_estimates.gibbs_read_count_samples.front()));
    } 
}


MinimumPathAbundanceEstimator::MinimumPathAbundanceEstimator(const uint32_t max_em_its, const double max_rel_em_conv, const uint32_t num_gibbs_samples, const uint32_t gibbs_thin_its, const double prob_precision) : PathAbundanceEstimator(max_em_its, max_rel_em_conv, num_gibbs_samples, gibbs_thin_its, prob_precision) {}

void MinimumPathAbundanceEstimator::estimate(PathClusterEstimates * path_cluster_estimates, const vector<ReadPathProbabilities> & cluster_probs, mt19937 * mt_rng) {

    if (!cluster_probs.empty()) {

        Utils::ColMatrixXd read_path_probs;
        Utils::ColVectorXd noise_probs;
        Utils::RowVectorXd read_counts;

        constructProbabilityMatrix(&read_path_probs, &noise_probs, &read_counts, cluster_probs, path_cluster_estimates->paths.size());      

        Utils::ColMatrixXb read_path_cover = Utils::ColMatrixXb::Zero(read_path_probs.rows(), read_path_probs.cols());
        Utils::RowVectorXd path_weights = Utils::RowVectorXd::Zero(read_path_probs.cols());

        for (size_t i = 0; i < read_path_probs.rows(); ++i) {

            if (Utils::doubleCompare(noise_probs(i), 1)) {

                read_counts(i) = 0;
            }

            for (auto & path_probs: cluster_probs.at(i).pathProbs()) {

                for (auto & path: path_probs.second) {

                    assert(path_probs.first > 0);

                    read_path_cover(i, path) = true;
                    path_weights(path) += log(path_probs.first) * read_counts(i);  
                }               
            }
        }

        path_weights *= -1;
        vector<uint32_t> min_path_cover = weightedMinimumPathCover(read_path_cover, read_counts, path_weights);

        if (!min_path_cover.empty()) {

            Utils::ColMatrixXd min_path_read_path_probs;
            Utils::ColVectorXd min_path_noise_probs;
            Utils::RowVectorXd min_path_read_counts;

            constructPartialProbabilityMatrix(&min_path_read_path_probs, &min_path_noise_probs, &min_path_read_counts, cluster_probs, min_path_cover, path_cluster_estimates->paths.size(), true);
            detractNoiseAndNormalizeProbabilityMatrix(&min_path_read_path_probs, &min_path_noise_probs, &min_path_read_counts);

            if (min_path_read_path_probs.rows() == 0) {

                assert(min_path_noise_probs.rows() == 0);
                assert(min_path_read_counts.cols() == 0);

                path_cluster_estimates->initEstimates(path_cluster_estimates->paths.size(), 0, true);
                return;
            }

            assert(min_path_read_path_probs.cols() >= 1);
            readCollapseProbabilityMatrix(&min_path_read_path_probs, &min_path_read_counts);

            const double total_min_path_read_counts = min_path_read_counts.sum();
            assert(total_min_path_read_counts > 0);

            PathClusterEstimates min_path_cluster_estimates;
            min_path_cluster_estimates.initEstimates(min_path_read_path_probs.cols(), 0, false);

            EMAbundanceEstimator(&min_path_cluster_estimates, min_path_read_path_probs, min_path_read_counts, total_min_path_read_counts);
            assert(min_path_cluster_estimates.abundances.cols() == min_path_cover.size());            

            path_cluster_estimates->initEstimates(path_cluster_estimates->paths.size(), 0, true);

            if (num_gibbs_samples > 0) {

                vector<CountSamples> * gibbs_read_count_samples = &(min_path_cluster_estimates.gibbs_read_count_samples);
                gibbs_read_count_samples->emplace_back(CountSamples());

                gibbs_read_count_samples->back().path_ids = min_path_cover;                
                gibbs_read_count_samples->back().samples = vector<vector<double> >(min_path_cluster_estimates.abundances.cols(), vector<double>());

                gibbsReadCountSampler(&min_path_cluster_estimates, min_path_read_path_probs, min_path_read_counts, total_min_path_read_counts, abundance_gibbs_gamma, mt_rng);
            }

            min_path_cluster_estimates.abundances *= total_min_path_read_counts;
            updateEstimates(path_cluster_estimates, min_path_cluster_estimates, min_path_cover, 1);    

        } else {

            path_cluster_estimates->initEstimates(path_cluster_estimates->paths.size(), 0, true);
        }

    } else {

        path_cluster_estimates->initEstimates(path_cluster_estimates->paths.size(), 0, true);
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


NestedPathAbundanceEstimator::NestedPathAbundanceEstimator(const uint32_t group_size_in, const uint32_t num_subset_samples_in, const bool infer_collapsed_in, const bool use_group_post_gibbs_in, const uint32_t max_em_its, const double max_rel_em_conv, const uint32_t num_gibbs_samples, const uint32_t gibbs_thin_its, const double prob_precision) : group_size(group_size_in), num_subset_samples(num_subset_samples_in), infer_collapsed(infer_collapsed_in), use_group_post_gibbs(use_group_post_gibbs_in), PathAbundanceEstimator(max_em_its, max_rel_em_conv, num_gibbs_samples, gibbs_thin_its, prob_precision) {}

void NestedPathAbundanceEstimator::estimate(PathClusterEstimates * path_cluster_estimates, const vector<ReadPathProbabilities> & cluster_probs, mt19937 * mt_rng) {

    if (infer_collapsed) {

        inferAbundancesCollapsedGroups(path_cluster_estimates, cluster_probs, mt_rng);
    
    } else {

        inferAbundancesIndependentGroups(path_cluster_estimates, cluster_probs, mt_rng);
    }
}

void NestedPathAbundanceEstimator::inferAbundancesIndependentGroups(PathClusterEstimates * path_cluster_estimates, const vector<ReadPathProbabilities> & cluster_probs, mt19937 * mt_rng) const {

    if (!cluster_probs.empty()) {

        auto path_groups = findPathGroups(path_cluster_estimates->paths);

        vector<vector<uint32_t> > path_subset_samples(num_subset_samples);

        for (auto & path_subset: path_subset_samples) {

            path_subset.reserve(path_groups.size() * group_size);
        }

        for (auto & group: path_groups) {    

            Utils::ColMatrixXd group_read_path_probs;
            Utils::ColVectorXd group_noise_probs;
            Utils::RowVectorXd group_read_counts;        

            constructPartialProbabilityMatrix(&group_read_path_probs, &group_noise_probs, &group_read_counts, cluster_probs, group, path_cluster_estimates->paths.size(), false);

            addNoiseAndNormalizeProbabilityMatrix(&group_read_path_probs, group_noise_probs);
            readCollapseProbabilityMatrix(&group_read_path_probs, &group_read_counts);

            group_noise_probs = group_read_path_probs.col(group_read_path_probs.cols() - 1);
            group_read_path_probs.conservativeResize(group_read_path_probs.rows(), group_read_path_probs.cols() - 1);

            vector<uint32_t> group_path_counts;
            group_path_counts.reserve(group.size());

            for (size_t i = 0; i < group.size(); ++i) {

                group_path_counts.emplace_back(path_cluster_estimates->paths.at(group.at(i)).source_count);
            }

            PathClusterEstimates group_path_cluster_estimates;

            if (use_group_post_gibbs) {

                estimatePathGroupPosteriorsGibbs(&group_path_cluster_estimates, group_read_path_probs, group_noise_probs, group_read_counts, group_path_counts, group_size, mt_rng);

            } else {

                if (group_size == 2) {

                    const double min_rel_likelihood = 1 / static_cast<double>(min_rel_likelihood_scaling * num_subset_samples); 
                    calculatePathGroupPosteriorsBounded(&group_path_cluster_estimates, group_read_path_probs, group_noise_probs, group_read_counts, group_path_counts, group_size, min_rel_likelihood);

                } else {

                    calculatePathGroupPosteriorsFull(&group_path_cluster_estimates, group_read_path_probs, group_noise_probs, group_read_counts, group_path_counts, group_size);                    
                }
            }

            sampleGroupPathIndices(&path_subset_samples, group_path_cluster_estimates, group, mt_rng);
        }

        spp::sparse_hash_map<vector<uint32_t>, uint32_t> clustered_path_subset_samples;

        for (auto & path_subset: path_subset_samples) {

            sort(path_subset.begin(), path_subset.end());

            auto clustered_path_subset_samples_it = clustered_path_subset_samples.emplace(path_subset, 0);
            clustered_path_subset_samples_it.first->second++;
        }

        inferPathSubsetAbundance(path_cluster_estimates, cluster_probs, mt_rng, clustered_path_subset_samples);

    } else {

        path_cluster_estimates->initEstimates(path_cluster_estimates->paths.size(), 0, true);
    }
}

void NestedPathAbundanceEstimator::inferAbundancesCollapsedGroups(PathClusterEstimates * path_cluster_estimates, const vector<ReadPathProbabilities> & cluster_probs, mt19937 * mt_rng) {

    if (!cluster_probs.empty()) {

        bool debug = false;     
        string debug_path = "";        

        for (auto & path: path_cluster_estimates->paths) {        

            if (
                path.name == "ENST00000393820.2_47" || 
                path.name == "ENST00000322723.8_227" || 
                path.name == "ENST00000322723.8_228" || 
                path.name == "ENST00000330899.4_62" || 
                path.name == "ENST00000360004.5_2543" || 
                path.name == "ENST00000360004.5_3301" || 
                path.name == "ENST00000300026.3_28" || 
                path.name == "ENST00000300026.3_53" || 
                path.name == "ENST00000374975.3_1122" || 
                path.name == "ENST00000374975.3_1517" || 
                path.name == "ENST00000557016.5_30" || 
                path.name == "ENST00000330459.7_80" || 
                path.name == "ENST00000234313.7_64" || 
                path.name == "ENST00000234313.7_110" || 
                path.name == "ENST00000366560.3_36"    
            ) {       

                debug = true;     
                debug_path = path.name;        

                break;        
            }      
        }

        auto path_source_groups = findPathSourceGroups(path_cluster_estimates->paths);

        if (debug) {        

            debug_mutex.lock();     

            cerr << "### " << debug_path << endl;      
            cerr << path_source_groups.first.size() << endl;      
            cerr << findPathGroups(path_cluster_estimates->paths).size() << endl;      

            uint32_t a = 0;       

            for (size_t i = 0; i < path_source_groups.first.size(); ++i) {        

                cerr << "#" <<  i << ": ";        

                for (auto & p: path_source_groups.first.at(i)) {      

                    cerr << path_cluster_estimates->paths.at(p).name << " ";      
                }      

                cerr << "\n" << path_source_groups.second.at(i) << endl;      
                a += path_source_groups.second.at(i);      
            }      

            cerr << a << endl;        
            cerr << cluster_probs.size() << " " << path_cluster_estimates->paths.size() << endl;       
        }

        Utils::ColMatrixXd group_read_path_probs;
        Utils::ColVectorXd group_noise_probs;
        Utils::RowVectorXd group_read_counts;  

        constructGroupedProbabilityMatrix(&group_read_path_probs, &group_noise_probs, &group_read_counts, cluster_probs, path_source_groups.first, path_cluster_estimates->paths.size());
        
        addNoiseAndNormalizeProbabilityMatrix(&group_read_path_probs, group_noise_probs);

        if (debug) {                

            cerr << "### " << debug_path << endl;
            cerr << group_read_path_probs.rows() << " " << group_read_path_probs.cols() << endl;      
        }

        readCollapseProbabilityMatrix(&group_read_path_probs, &group_read_counts);

        if (debug) {                

            cerr << "### " << debug_path << endl;      
            cerr << group_read_path_probs.rows() << " " << group_read_path_probs.cols() << endl;           
         }

        group_noise_probs = group_read_path_probs.col(group_read_path_probs.cols() - 1);
        group_read_path_probs.conservativeResize(group_read_path_probs.rows(), group_read_path_probs.cols() - 1);

        PathClusterEstimates group_path_cluster_estimates;

        if (use_group_post_gibbs) {

            estimatePathGroupPosteriorsGibbs(&group_path_cluster_estimates, group_read_path_probs, group_noise_probs, group_read_counts, path_source_groups.second, group_size, mt_rng);

        } else {

            if (group_size == 2) {

                const double min_rel_likelihood = 1 / static_cast<double>(min_rel_likelihood_scaling * num_subset_samples); 
                calculatePathGroupPosteriorsBounded(&group_path_cluster_estimates, group_read_path_probs, group_noise_probs, group_read_counts, path_source_groups.second, group_size, min_rel_likelihood);

            } else {

                calculatePathGroupPosteriorsFull(&group_path_cluster_estimates, group_read_path_probs, group_noise_probs, group_read_counts, path_source_groups.second, group_size);                    
            }
        }

        if (debug) {                

            cerr << "### " << debug_path << endl;      
            cerr << group_path_cluster_estimates.posteriors.size() << endl;       
            cerr << group_path_cluster_estimates.path_group_sets.size() << endl;       

            for (size_t i = 0; i < group_path_cluster_estimates.path_group_sets.size(); ++i) {        

                if (group_path_cluster_estimates.posteriors.at(i) > pow(10, -4)) {        

                    cerr << group_path_cluster_estimates.posteriors.at(i) << ": " << group_path_cluster_estimates.path_group_sets.at(i) << endl;      
                }      
            }               

        }

        spp::sparse_hash_map<vector<uint32_t>, uint32_t> path_subset_samples;
        samplePathSubsetIndices(&path_subset_samples, group_path_cluster_estimates, path_source_groups.first, mt_rng);

        if (debug) {                

            cerr << "### " << debug_path << endl;
            cerr << path_subset_samples.size() << endl;       

            for (auto & bla: path_subset_samples) {       

                for (auto & bla2: bla.first) {   

                    cerr << path_cluster_estimates->paths.at(bla2).name << "," << bla2 << " ";    
                }

                cerr << ": " << bla.second << endl;      
            }      

            debug_mutex.unlock();     
         }

        inferPathSubsetAbundance(path_cluster_estimates, cluster_probs, mt_rng, path_subset_samples);

    } else {

        path_cluster_estimates->initEstimates(path_cluster_estimates->paths.size(), 0, true);
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

pair<vector<vector<uint32_t> >, vector<uint32_t> > NestedPathAbundanceEstimator::findPathSourceGroups(const vector<PathInfo> & paths) const {

    spp::sparse_hash_map<uint32_t, vector<uint32_t> > source_id_paths;

    for (size_t i = 0; i < paths.size(); ++i) {

        for (auto & id: paths.at(i).source_ids) {

            auto source_id_paths_it = source_id_paths.emplace(id, vector<uint32_t>());
            source_id_paths_it.first->second.emplace_back(i);
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

void NestedPathAbundanceEstimator::sampleGroupPathIndices(vector<vector<uint32_t> > * path_subset_samples, const PathClusterEstimates & group_path_cluster_estimates, const vector<uint32_t> & group, mt19937 * mt_rng) const {

    assert(group_path_cluster_estimates.posteriors.size() == group_path_cluster_estimates.path_group_sets.size());
    discrete_distribution<uint32_t> path_group_set_sampler(group_path_cluster_estimates.posteriors.begin(), group_path_cluster_estimates.posteriors.end());

    for (size_t i = 0; i < num_subset_samples; ++i) {

        vector<uint32_t> path_group_set = group_path_cluster_estimates.path_group_sets.at(path_group_set_sampler(*mt_rng));

        assert(!path_group_set.empty());
        assert(path_group_set.size() == group_size);

        sort(path_group_set.begin(), path_group_set.end());

        for (auto & path_group: path_group_set) {

            path_subset_samples->at(i).emplace_back(group.at(path_group));
        }
    }
}

void NestedPathAbundanceEstimator::samplePathSubsetIndices(spp::sparse_hash_map<vector<uint32_t>, uint32_t> * path_subset_samples, const PathClusterEstimates & group_path_cluster_estimates, const vector<vector<uint32_t> > & path_groups, mt19937 * mt_rng) const {

    assert(group_path_cluster_estimates.posteriors.size() == group_path_cluster_estimates.path_group_sets.size());
    discrete_distribution<uint32_t> path_group_set_sampler(group_path_cluster_estimates.posteriors.begin(), group_path_cluster_estimates.posteriors.end());

    vector<uint32_t> path_group_set_sample_counts(group_path_cluster_estimates.path_group_sets.size(), 0);

    for (size_t i = 0; i < num_subset_samples; ++i) {

        path_group_set_sample_counts.at(path_group_set_sampler(*mt_rng)) += 1;
    }

    for (size_t i = 0; i < path_group_set_sample_counts.size(); ++i) {

        if (path_group_set_sample_counts.at(i) > 0) {

            auto path_group_set = group_path_cluster_estimates.path_group_sets.at(i);

            assert(!path_group_set.empty());
            assert(path_group_set.size() == group_size);

            vector<uint32_t> path_subset;

            for (auto & group: path_group_set) {

                for (auto & path: path_groups.at(group)) {

                    path_subset.emplace_back(path);
                }
            }

            sort(path_subset.begin(), path_subset.end());
            auto path_subset_samples_it = path_subset_samples->emplace(path_subset, 0);

            path_subset_samples_it.first->second += path_group_set_sample_counts.at(i);
        }
    }
}

void NestedPathAbundanceEstimator::inferPathSubsetAbundance(PathClusterEstimates * path_cluster_estimates, const vector<ReadPathProbabilities> & cluster_probs, mt19937 * mt_rng, const spp::sparse_hash_map<vector<uint32_t>, uint32_t> & path_subset_samples) const {

    path_cluster_estimates->initEstimates(path_cluster_estimates->paths.size(), 0, true);

    spp::sparse_hash_map<vector<uint32_t>, uint32_t> subset_path_group_samples;    

    for (auto & path_subset: path_subset_samples) {

        assert(!path_subset.first.empty());
        assert(path_subset.second > 0);

        spp::sparse_hash_map<uint32_t, vector<uint32_t> > subset_path_group;
        
        for (auto & path: path_subset.first) {

            auto subset_path_group_it = subset_path_group.emplace(path_cluster_estimates->paths.at(path).group_id, vector<uint32_t>());
            subset_path_group_it.first->second.emplace_back(path);
        }

        for (auto & path_group: subset_path_group) {

            auto subset_path_group_samples_it = subset_path_group_samples.emplace(path_group.second, 0);
            subset_path_group_samples_it.first->second += path_subset.second;
        }

        vector<uint32_t> collapsed_path_subset;
        collapsed_path_subset.reserve(path_subset.first.size());

        collapsed_path_subset.emplace_back(path_subset.first.front());

        for (size_t i = 1; i < path_subset.first.size(); ++i) {

            if (path_subset.first.at(i) != collapsed_path_subset.back()) {

                collapsed_path_subset.emplace_back(path_subset.first.at(i));
            }
        }

        Utils::ColMatrixXd subset_read_path_probs;
        Utils::ColVectorXd subset_noise_probs;
        Utils::RowVectorXd subset_read_counts;

        constructPartialProbabilityMatrix(&subset_read_path_probs, &subset_noise_probs, &subset_read_counts, cluster_probs, collapsed_path_subset, path_cluster_estimates->paths.size(), true);
        detractNoiseAndNormalizeProbabilityMatrix(&subset_read_path_probs, &subset_noise_probs, &subset_read_counts);

        if (subset_read_path_probs.rows() == 0) {

            assert(subset_noise_probs.rows() == 0);
            assert(subset_read_counts.cols() == 0);

            path_cluster_estimates->initEstimates(path_cluster_estimates->paths.size(), 0, true);
            return;
        }

        assert(subset_read_path_probs.cols() >= 1);
        readCollapseProbabilityMatrix(&subset_read_path_probs, &subset_read_counts);

        const double total_subset_read_counts = subset_read_counts.sum();
        assert(total_subset_read_counts > 0);

        PathClusterEstimates subset_path_cluster_estimates;
        subset_path_cluster_estimates.initEstimates(subset_read_path_probs.cols(), 0, false);

        EMAbundanceEstimator(&subset_path_cluster_estimates, subset_read_path_probs, subset_read_counts, total_subset_read_counts);
        assert(subset_path_cluster_estimates.abundances.cols() == collapsed_path_subset.size());            

        if (num_gibbs_samples > 0) {

            vector<CountSamples> * gibbs_read_count_samples = &(subset_path_cluster_estimates.gibbs_read_count_samples);
            gibbs_read_count_samples->emplace_back(CountSamples());

            gibbs_read_count_samples->back().path_ids = collapsed_path_subset;            
            gibbs_read_count_samples->back().samples = vector<vector<double> >(subset_path_cluster_estimates.abundances.cols(), vector<double>());

            for (uint32_t i = 0; i < path_subset.second; ++i) {

                gibbsReadCountSampler(&subset_path_cluster_estimates, subset_read_path_probs, subset_read_counts, total_subset_read_counts, abundance_gibbs_gamma, mt_rng);
            }
        }

        subset_path_cluster_estimates.abundances *= total_subset_read_counts;
        updateEstimates(path_cluster_estimates, subset_path_cluster_estimates, collapsed_path_subset, path_subset.second);
    }

    assert(path_cluster_estimates->posteriors.empty());
    assert(path_cluster_estimates->path_group_sets.empty());

    path_cluster_estimates->posteriors.reserve(subset_path_group_samples.size());
    path_cluster_estimates->path_group_sets.reserve(subset_path_group_samples.size());

    for (auto & path_group_sample: subset_path_group_samples) {

        assert(path_group_sample.first.size() <= group_size);

        path_cluster_estimates->posteriors.emplace_back(path_group_sample.second / static_cast<double>(num_subset_samples));
        path_cluster_estimates->path_group_sets.emplace_back(path_group_sample.first);
    }    

    for (size_t i = 0; i < path_cluster_estimates->abundances.cols(); ++i) {

        path_cluster_estimates->abundances(0, i) /= static_cast<double>(num_subset_samples);
    }
}

