
#ifndef FERSKEN_SRC_PATHABUNDANCEESTIMATOR_HPP
#define FERSKEN_SRC_PATHABUNDANCEESTIMATOR_HPP

#include <vector>
#include <random>

#include <Eigen/Dense>

#include "path_abundances.hpp"
#include "read_path_probabilities.hpp"
#include "utils.hpp"

using namespace std;


namespace Eigen {

    typedef Eigen::Matrix<double, Eigen::Dynamic, 1, Eigen::ColMajor> ColVectorXd;
    
    typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> RowMatrixXd;
    typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor> ColMatrixXd;
}

class PathAbundanceEstimator {

    public:

        PathAbundanceEstimator(const double min_abundance_in, const uint32_t rng_seed);
        virtual ~PathAbundanceEstimator() {};

        virtual Abundances inferPathClusterAbundances(const vector<pair<ReadPathProbabilities, uint32_t> > & cluster_probs, const uint32_t num_paths) const = 0;

    protected: 

        const double min_abundance;
        mt19937 mt_rng;

        void nullifyAbundances(Abundances * abundances) const;
        void removeNoiseAndRenormalizeAbundances(Abundances * abundances) const;
};

class SimplePathAbundanceEstimator : public PathAbundanceEstimator {

    public:

        SimplePathAbundanceEstimator(const uint32_t max_em_iterations_in, const double min_abundance, const uint32_t rng_seed = 0);
        ~SimplePathAbundanceEstimator() {};

        Abundances inferPathClusterAbundances(const vector<pair<ReadPathProbabilities, uint32_t> > & cluster_probs, const uint32_t num_paths) const;

    protected:

        const uint32_t max_em_iterations;
    
        void expectationMaximizationEstimator(Abundances * abundances, const Eigen::RowMatrixXd & read_path_probs, const Eigen::RowVectorXi & read_counts) const;
};

class DiploidPathAbundanceEstimator : public SimplePathAbundanceEstimator {

    public:

        DiploidPathAbundanceEstimator(const uint32_t max_diploid_iterations_in, const uint32_t max_em_iterations, const double min_abundance, const uint32_t rng_seed);
        ~DiploidPathAbundanceEstimator() {};

        Abundances inferPathClusterAbundances(const vector<pair<ReadPathProbabilities, uint32_t> > & cluster_probs, const uint32_t num_paths) const;

    private:
    
        const uint32_t max_diploid_iterations;
};


#endif
