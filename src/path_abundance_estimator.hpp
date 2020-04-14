
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

        virtual Abundances inferPathClusterAbundances(const vector<pair<ReadPathProbabilities, uint32_t> > & cluster_probs, const uint32_t num_paths);

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

        Abundances inferPathClusterAbundances(const vector<pair<ReadPathProbabilities, uint32_t> > & cluster_probs, const uint32_t num_paths);

        vector<uint32_t> weightedMinimumPathCover(const Eigen::ColMatrixXb & read_path_cover, const Eigen::RowVectorXui & read_counts, const Eigen::RowVectorXd & path_weights);
};


#endif
