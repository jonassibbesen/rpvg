
#ifndef RPVG_SRC_PATHABUNDANCEESTIMATOR_HPP
#define RPVG_SRC_PATHABUNDANCEESTIMATOR_HPP

#include <vector>
#include <random>
#include <mutex>

#include <Eigen/Dense>

#include "path_abundances.hpp"
#include "read_path_probabilities.hpp"
#include "utils.hpp"

using namespace std;


class PathAbundanceEstimator {

    public:

        PathAbundanceEstimator(const uint32_t max_em_its_in, const double min_read_count_in);
        virtual ~PathAbundanceEstimator() {};

        virtual PathAbundances inferPathClusterAbundances(const vector<pair<ReadPathProbabilities, uint32_t> > & cluster_probs, const vector<Path> & cluster_paths);

    protected: 

        const uint32_t max_em_its;
        const double min_read_count;

        void expectationMaximizationEstimator(Abundances * abundances, const Eigen::ColMatrixXd & read_path_probs, const Eigen::RowVectorXui & read_counts) const;
        void removeNoiseAndRenormalizeAbundances(Abundances * abundances) const;    
};

class MinimumPathAbundanceEstimator : public PathAbundanceEstimator {

    public:

        MinimumPathAbundanceEstimator(const uint32_t max_em_its, const double min_read_count);
        ~MinimumPathAbundanceEstimator() {};

        PathAbundances inferPathClusterAbundances(const vector<pair<ReadPathProbabilities, uint32_t> > & cluster_probs, const vector<Path> & cluster_paths);

        vector<uint32_t> weightedMinimumPathCover(const Eigen::ColMatrixXb & read_path_cover, const Eigen::RowVectorXui & read_counts, const Eigen::RowVectorXd & path_weights);
};

class NestedPathAbundanceEstimator : public PathAbundanceEstimator {

    public:

        NestedPathAbundanceEstimator(const uint32_t num_nested_its_in, const uint32_t ploidy_in, const uint32_t rng_seed, const uint32_t max_em_its, const double min_read_count);
        ~NestedPathAbundanceEstimator() {};

        PathAbundances inferPathClusterAbundances(const vector<pair<ReadPathProbabilities, uint32_t> > & cluster_probs, const vector<Path> & cluster_paths);

    private:

        const uint32_t num_nested_its;
        const uint32_t ploidy;

        mt19937 mt_rng;

        spp::sparse_hash_map<string, uint32_t> grouping;
};

namespace std {

    template<> 
    struct hash<vector<uint32_t> >
    {
        size_t operator()(vector<uint32_t> const & vec) const
        {
            size_t seed = 0;

            for (auto & val: vec) {

                spp::hash_combine(seed, val);
            }

            return seed;
        }
    };
}

 
#endif
