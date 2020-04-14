
#ifndef FERSKEN_SRC_PATHABUNDANCEESTIMATOR_HPP
#define FERSKEN_SRC_PATHABUNDANCEESTIMATOR_HPP

#include <vector>
#include <random>

#include <Eigen/Dense>

#include "path_abundances.hpp"
#include "read_path_probabilities.hpp"
#include "utils.hpp"

using namespace std;


class PathAbundanceEstimator {

    public:

        PathAbundanceEstimator(const uint32_t max_em_its_in, const double min_abundance_in);
        virtual ~PathAbundanceEstimator() {};

        virtual PathAbundances inferPathClusterAbundances(const vector<pair<ReadPathProbabilities, uint32_t> > & cluster_probs, const vector<Path> & cluster_paths);

    protected: 

        const uint32_t max_em_its;
        const double min_abundance;

        void expectationMaximizationEstimator(Abundances * abundances, const Eigen::ColMatrixXd & read_path_probs, const Eigen::RowVectorXui & read_counts) const;
        void removeNoiseAndRenormalizeAbundances(Abundances * abundances) const;    
};

class MinimumPathAbundanceEstimator : public PathAbundanceEstimator {

    public:

        MinimumPathAbundanceEstimator(const uint32_t max_em_its, const double min_abundance);
        ~MinimumPathAbundanceEstimator() {};

        PathAbundances inferPathClusterAbundances(const vector<pair<ReadPathProbabilities, uint32_t> > & cluster_probs, const vector<Path> & cluster_paths);

        vector<uint32_t> weightedMinimumPathCover(const Eigen::ColMatrixXb & read_path_cover, const Eigen::RowVectorXui & read_counts, const Eigen::RowVectorXd & path_weights);
};

class GroupedPathAbundanceEstimator : public PathAbundanceEstimator {

    public:

        GroupedPathAbundanceEstimator(const uint32_t num_group_its_in, const uint32_t ploidy_in, const uint32_t rng_seed, const uint32_t max_em_its, const double min_abundance);
        ~GroupedPathAbundanceEstimator() {};

        PathAbundances inferPathClusterAbundances(const vector<pair<ReadPathProbabilities, uint32_t> > & cluster_probs, const vector<Path> & cluster_paths);

    private:

        const uint32_t num_group_its;
        const uint32_t ploidy;

        mt19937 mt_rng;

        spp::sparse_hash_map<string, uint32_t> grouping;
};

 
#endif
