
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

    typedef Eigen::Matrix<uint32_t, Eigen::Dynamic, 1, Eigen::ColMajor> ColVectorXui;
    typedef Eigen::Matrix<double, Eigen::Dynamic, 1, Eigen::ColMajor> ColVectorXd;

    typedef Eigen::Matrix<uint32_t, 1, Eigen::Dynamic, Eigen::RowMajor> RowVectorXui;
    typedef Eigen::Matrix<double, 1, Eigen::Dynamic, Eigen::RowMajor> RowVectorXd;
    
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
    
        void expectationMaximizationEstimator(Abundances * abundances, const Eigen::ColMatrixXd & read_path_probs, const Eigen::RowVectorXui & read_counts) const;
};

class MinimumPathAbundanceEstimator : public SimplePathAbundanceEstimator {

    public:

        MinimumPathAbundanceEstimator(const uint32_t num_path_iterations_in, const uint32_t max_em_iterations, const double min_abundance, const uint32_t rng_seed);
        ~MinimumPathAbundanceEstimator() {};

        Abundances inferPathClusterAbundances(const vector<pair<ReadPathProbabilities, uint32_t> > & cluster_probs, const uint32_t num_paths);

    private:
    
        const uint32_t num_path_iterations;

        vector<uint32_t> sampleMinimumPathCover(const Eigen::ColMatrixXd & read_path_noise_log_probs, const Eigen::RowVectorXui & read_counts);
};

class DiploidPathAbundanceEstimator : public SimplePathAbundanceEstimator {

    public:

        DiploidPathAbundanceEstimator(const uint32_t num_diploid_iterations_in, const uint32_t max_em_iterations, const double min_abundance, const uint32_t rng_seed);
        ~DiploidPathAbundanceEstimator() {};

        Abundances inferPathClusterAbundances(const vector<pair<ReadPathProbabilities, uint32_t> > & cluster_probs, const uint32_t num_paths) const;

    private:
    
        const uint32_t num_diploid_iterations;
};


#endif
